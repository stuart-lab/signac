"""Helpers for querying the AlphaGenome Atlas API from Signac (via reticulate).

The Atlas API can respond with UNAVAILABLE when a key exceeds its request
rate. The functions here wrap the official alphagenome Atlas client with
per-request retries (exponential backoff with jitter, shared across worker
threads) so that a region can be retrieved despite intermittent throttling.
"""

import concurrent.futures
import random
import threading
import time

import anndata
import grpc
import numpy as np
import pandas as pd

from alphagenome.data import genome

# The Atlas interval endpoint rejects requests wider than 32 bp.
CHUNK_SIZE = 32
_RETRYABLE = ("UNAVAILABLE", "RESOURCE_EXHAUSTED", "DEADLINE_EXCEEDED")


class _Throttle:
  """Shared backoff state so that all workers pause after a throttled call."""

  def __init__(self):
    self._lock = threading.Lock()
    self._resume_at = 0.0

  def wait(self):
    delay = self._resume_at - time.monotonic()
    if delay > 0:
      time.sleep(delay)

  def backoff(self, attempt, max_delay=60.0):
    delay = min(max_delay, (2.0**attempt)) * (0.5 + random.random())
    with self._lock:
      self._resume_at = max(self._resume_at, time.monotonic() + delay)
    return delay


def _is_retryable(error):
  if isinstance(error, TimeoutError):
    return True
  cause = error
  while cause is not None:
    if isinstance(cause, grpc.RpcError) and callable(getattr(cause, "code", None)):
      return cause.code() in (
          grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.RESOURCE_EXHAUSTED,
          grpc.StatusCode.DEADLINE_EXCEEDED,
      )
    cause = cause.__cause__
  text = str(error)
  return any(code in text for code in _RETRYABLE)


def _with_retries(fn, throttle, max_retries):
  """Calls fn(), retrying throttled/unavailable responses with backoff."""
  attempt = 0
  while True:
    throttle.wait()
    try:
      return fn()
    except ValueError:
      # invalid argument / not found: not retryable
      raise
    except Exception as error:  # pylint: disable=broad-except
      if not _is_retryable(error) or attempt >= max_retries:
        raise
      throttle.backoff(attempt)
      attempt += 1


def _concat(results):
  """Concatenates a list of AnnData objects along the observation axis."""
  results = [r for r in results if r is not None and r.n_obs > 0]
  if not results:
    return None
  if len(results) == 1:
    return results[0]
  offset = 0
  for res in results:
    res.obs.index = [str(offset + i) for i in range(res.n_obs)]
    offset += res.n_obs
  combined = anndata.concat(results, join="outer")
  combined.var = results[0].var
  return combined


def _list_scores_direct(client, interval, scorer, filter_query, field_mask):
  """Calls ListDenseVariantScores for one chunk, following pagination."""
  from alphagenome.atlas import atlas as atlas_module
  from alphagenome.protos import atlas_service_pb2, dna_model_pb2

  request = atlas_service_pb2.ListDenseVariantScoresRequest(
      interval=interval.to_proto(),
      organism=dna_model_pb2.ORGANISM_HOMO_SAPIENS,
      filter=filter_query,
  )
  metadata = [*client._metadata, ("x-goog-fieldmask", field_mask)]  # pylint: disable=protected-access
  scores = []
  while True:
    with atlas_module.handle_rpc_error():
      response = client._stub.ListDenseVariantScores(  # pylint: disable=protected-access
          request, metadata=metadata
      )
    scores.extend(response.variant_scores)
    if not response.next_page_token:
      return scores
    request.page_token = response.next_page_token


def query_interval_chunked(
    client,
    chromosome,
    start,
    end,
    scorer,
    max_workers=4,
    max_retries=5,
):
  """Queries the Atlas interval endpoint in 32 bp chunks with retries.

  When given a real AtlasClient the requests are issued directly against the
  gRPC stub, so the scorer metadata is fetched once for the whole region
  rather than once per chunk (as AtlasClient.query_interval does), and the
  scores of all chunks are converted to a single AnnData at the end. Clients
  without a stub (e.g. test doubles) are queried through query_interval.

  Args:
    client: An alphagenome.atlas.atlas.AtlasClient.
    chromosome: Chromosome name (e.g. 'chr1').
    start: 1-based start position (inclusive).
    end: 1-based end position (inclusive).
    scorer: Name of the variant scorer to request.
    max_workers: Number of chunks requested concurrently.
    max_retries: Maximum number of retries per chunk after throttling.

  Returns:
    An AnnData with one row per variant, or None if nothing was returned.
  """
  start = int(start)
  end = int(end)
  throttle = _Throttle()
  chunks = []
  pos = start - 1  # 0-based half-open
  while pos < end:
    chunks.append(genome.Interval(chromosome, pos, min(pos + CHUNK_SIZE, end)))
    pos += CHUNK_SIZE

  direct = hasattr(client, "_stub") and hasattr(client, "_metadata")
  if direct:
    from alphagenome.atlas import atlas as atlas_module
    from alphagenome.atlas import atlas_utils

    filter_query = atlas_utils.build_filter(requested_scorers=[scorer])
    field_mask = ",".join(atlas_module._LIST_DENSE_VARIANT_SCORES_FIELD_MASKS)  # pylint: disable=protected-access
    scorer_metadata = _with_retries(
        lambda: atlas_module._filter_scorer_metadata(client.scorer_metadata()),  # pylint: disable=protected-access
        throttle,
        max_retries,
    )

    def _fetch(interval):
      return _with_retries(
          lambda: _list_scores_direct(
              client, interval, scorer, filter_query, field_mask
          ),
          throttle,
          max_retries,
      )

    scores = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=int(max_workers)
    ) as executor:
      for res in executor.map(_fetch, chunks):
        scores.extend(res)
    if not scores:
      return None
    results = atlas_module.convert_variant_scores_to_anndata(
        scores, scorer_metadata
    )
    return results.get(scorer)

  def _fetch(interval):
    def _call():
      res = client.query_interval(
          interval,
          requested_scorers=[scorer],
          progress_bar=False,
          max_workers=1,
      )
      return res.get(scorer)

    return _with_retries(_call, throttle, max_retries)

  results = []
  with concurrent.futures.ThreadPoolExecutor(
      max_workers=int(max_workers)
  ) as executor:
    for res in executor.map(_fetch, chunks):
      results.append(res)
  return _concat(results)


def variants_to_frame(ad):
  """Extracts the variants of an Atlas AnnData result into a DataFrame.

  Args:
    ad: An AnnData returned by the Atlas client, with genome.Variant objects
      in obs['variant'].

  Returns:
    A DataFrame with columns chromosome, position (1-based), ref and alt, in
    the same order as the rows of ad.
  """
  variants = list(ad.obs["variant"])
  return pd.DataFrame({
      "chromosome": [v.chromosome for v in variants],
      "position": np.asarray([v.position for v in variants], dtype=np.int64),
      "ref": [v.reference_bases for v in variants],
      "alt": [v.alternate_bases for v in variants],
  })
