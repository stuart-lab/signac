// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <zlib.h>
#include <string.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP groupCommand(std::string fragments) {
  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(fragments.c_str(), "rb");

  // return empty list if it can't find the file
  if (fileHandler == NULL) {
    Rcpp::Rcerr << "can't open file" << std::flush;
    return (Rcpp::DataFrame::create());
  }

  std::unordered_map<std::string, size_t> group_hash;

  char* cb_char;
  size_t column_idx {3};
  size_t line_counter {0};
  uint32_t buffer_length = 256;
  char *buffer = new char[buffer_length];
  while(gzgets(fileHandler, buffer, buffer_length) !=0 ){
    cb_char = strtok ( buffer, "\t" );
    for (auto i=1; i<=column_idx; i++) {
      cb_char = strtok (NULL, "\t");
    }

    std::string cb(cb_char);
    group_hash[cb] += 1;

    line_counter += 1;
    if (line_counter % 10000000 == 0) {
      Rcpp::Rcerr << "\r                                                  "
                  << std::flush;
    }
    if (line_counter % 1000000 == 0) {
      Rcpp::Rcerr << "\rDone Processing " << line_counter / 1000000
                  << " million lines";
    }


  }

  //Cleanup
  gzclose(fileHandler);

  std::vector<std::string> keys;
  keys.reserve(group_hash.size());

  std::vector<size_t> vals;
  vals.reserve(group_hash.size());

  for(auto kv : group_hash) {
    keys.push_back(kv.first);
    vals.push_back(kv.second);
  }

  Rcpp::DataFrame cb_groups = Rcpp::DataFrame::create(Rcpp::Named("CB") = keys,
                                                      Rcpp::Named("count") = vals);

  return (cb_groups);
}
