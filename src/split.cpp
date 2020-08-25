#include <Rcpp.h>
#include <zlib.h>
#include <iostream>
#include <fstream>


// [[Rcpp::export]]
int splitFragments(
    std::string fragments,
    std::vector<std::string> cells,
    std::vector<std::string> idents,
    std::vector<std::string> unique_idents,
    std::string outdir,
    std::string suffix,
    int buffer_length,
    bool append = false,
    bool verbose = true
) {
  // cells and idents should be the same length and order
  if (cells.size() != idents.size()) {
    return 1;
  }

  // construct map
  // cell is the key, value is the ident
  std::map<std::string, std::string> cellmap;
  for (size_t i = 0; i < cells.size(); ++i) {
    cellmap[cells[i]] = idents[i];
  }

  // opening gzipped compressed stream
  gzFile ifileHandler = gzopen(fragments.c_str(), "rb");

  // open one output file for each unique ident
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < unique_idents.size(); i++) {
    std::string fileName = outdir + unique_idents[i] + suffix + ".bed";
    if (append) {
      streams.emplace_back(std::ofstream{ fileName, std::ios_base::app});
    } else {
      streams.emplace_back(std::ofstream{ fileName });
    }
  }

  // return 1 if it can't find the file
  if (ifileHandler == NULL) {
    Rcpp::Rcerr << "can't open file" << std::flush;
    return 1;
  }

  // C based buffered string parsing
  char* cb_char;
  size_t line_counter {1};
  char *buffer = new char[buffer_length];

  size_t num_whitelist_cells {0};
  {

    if (verbose) {
      num_whitelist_cells = cellmap.size();

      Rcpp::Rcerr << "Keeping " << num_whitelist_cells
                  << " cell barcodes"
                  << std::endl << std::flush;

      Rcpp::Rcerr << "Splitting into " << unique_idents.size()
                  << " files"
                  << std::endl << std::flush;
    }
  }

  // char * to string extraction
  std::string cb_seq, line_seq;
  cb_seq.reserve(32);
  line_seq.reserve(buffer_length);

  // looping over the fragments file
  while(gzgets(ifileHandler, buffer, buffer_length) !=0 ){
    line_seq.clear();
    line_seq.append(buffer);

    cb_char = strtok ( buffer, "\t" );

    for (auto i=1; i<=3; i++) {
      cb_char = strtok (NULL, "\t");

      if(i == 3) {
        cb_seq.clear();
        cb_seq.append(cb_char);
        auto it = cellmap.find(cb_seq);
        if (it != cellmap.end()) {
          // look up ident of matched cell
          std::string cellident = it->second;

          // find which file to write to
          std::vector<std::string>::iterator iter = std::find(unique_idents.begin(), unique_idents.end(), cellident);
          int index = std::distance(unique_idents.begin(), iter);

          // write to correct stream for given ident
          streams[index] << line_seq.c_str();
        }
      }
    }

    line_counter += 1;
    bool is_ten_mil = line_counter % 10000000 == 0;
    if (verbose) {
      if (is_ten_mil) {
        Rcpp::Rcerr << "\r                                                  ";
      }

      if (line_counter % 1000000 == 0) {
        Rcpp::Rcerr << "\rDone Processing " << line_counter / 1000000
                    << " million lines";
      }
    }

    if (is_ten_mil) {
      Rcpp::checkUserInterrupt();
    }
  }

  // Cleanup
  gzclose(ifileHandler);
  for (int i = 0; i < streams.size(); i++) {
    streams[i].close();
  }

  return 0;
}
