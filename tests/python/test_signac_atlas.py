"""Offline unit tests for the Signac Atlas helper."""
import io
import unittest
from unittest.mock import patch

import anndata
import grpc
import numpy as np
import pandas as pd
import signac_atlas as atlas


class RpcFailure(grpc.RpcError):
    def __init__(self, status):
        self.status = status

    def code(self):
        return self.status


class DirectStub:
    """Fake gRPC stub answering ListDenseVariantScores with two pages."""

    def __init__(self, scorer):
        self.scorer = scorer
        self.requests = []

    def ListDenseVariantScores(self, request, metadata=None):  # noqa: N802
        from alphagenome.data import genome
        from alphagenome.protos import atlas_service_pb2
        self.requests.append((request.interval.start, request.interval.end,
                              request.page_token, tuple(metadata or ())))
        interval = genome.Interval.from_proto(request.interval)
        positions = list(range(interval.start + 1, interval.end + 1))
        if request.page_token == '':
            page, token = positions[:1], 'more' if len(positions) > 1 else ''
        else:
            page, token = positions[1:], ''
        variant_scores = []
        for pos in page:
            for alt in 'CGT':
                variant = genome.Variant(interval.chromosome, pos, 'A', alt)
                score = atlas_service_pb2.DenseVariantScore(
                    variant_scorer=atlas_service_pb2.VariantScorerInfo(
                        name=self.scorer, is_signed=False),
                    shape=[1, 2],
                    scores=np.array([[pos, 0.5]], dtype=np.float32).tobytes())
                variant_scores.append(atlas_service_pb2.DenseVariantScores(
                    variant=variant.to_proto(), scores=[score]))
        return atlas_service_pb2.ListDenseVariantScoresResponse(
            variant_scores=variant_scores, next_page_token=token)


class DirectClient:
    """Fake AtlasClient exposing the private attributes the direct path uses."""

    def __init__(self, scorer='AVI_SCORE_FEATURE_IMPORTANCE'):
        self.scorer = scorer
        self._stub = DirectStub(scorer)
        self._metadata = [('x-goog-api-key', 'key')]
        self.metadata_calls = 0

    def scorer_metadata(self):
        from alphagenome.atlas.atlas import ScorerMetadata
        self.metadata_calls += 1
        return {self.scorer: ScorerMetadata(
            name=self.scorer, is_signed=False,
            track_metadata=pd.DataFrame({'name': ['F1', 'F2']}, index=['0', '1']))}


class AtlasRegressions(unittest.TestCase):
    def test_direct_stub_path(self):
        client = DirectClient()
        result = atlas.query_interval_chunked(
            client, 'chr1', 1001, 1040, client.scorer,
            max_workers=1, max_retries=0)
        # two chunks (32 + 8 bp), each fetched in two pages
        self.assertEqual([r[:3] for r in client._stub.requests],
                         [(1000, 1032, ''), (1000, 1032, 'more'),
                          (1032, 1040, ''), (1032, 1040, 'more')])
        self.assertTrue(all('x-goog-fieldmask' in dict(r[3]) for r in
                            client._stub.requests))
        self.assertEqual(client.metadata_calls, 1)
        self.assertEqual(result.n_obs, 40 * 3)
        self.assertEqual(list(result.var['name']), ['F1', 'F2'])
        positions = sorted(v.position for v in result.obs['variant'])
        self.assertEqual(positions, sorted(list(range(1001, 1041)) * 3))
        self.assertTrue(np.allclose(
            sorted(result.X[:, 0]), sorted(positions)))
        frame = atlas.variants_to_frame(result)
        self.assertEqual(set(frame['alt']), {'C', 'G', 'T'})
        self.assertTrue((frame['ref'] == 'A').all())

    def test_direct_stub_path_returns_none_for_missing_scorer(self):
        client = DirectClient()
        result = atlas.query_interval_chunked(
            client, 'chr1', 1001, 1002, 'OTHER', max_workers=1, max_retries=0)
        self.assertIsNone(result)

    def test_timeout_retry(self):
        calls = []

        def fetch():
            calls.append(1)
            if len(calls) == 1:
                raise TimeoutError('Deadline expired before operation could complete.')
            return 'ok'

        with patch.object(atlas._Throttle, 'backoff'):
            self.assertEqual(atlas._with_retries(fetch, atlas._Throttle(), 2), 'ok')
        self.assertEqual(len(calls), 2)
        self.assertTrue(atlas._is_retryable(RpcFailure(grpc.StatusCode.UNAVAILABLE)))
        self.assertFalse(atlas._is_retryable(RpcFailure(grpc.StatusCode.PERMISSION_DENIED)))


def run_suite():
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(AtlasRegressions))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return True


if __name__ == '__main__':
    run_suite()
