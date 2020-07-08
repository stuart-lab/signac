#include <Rcpp.h>
#include <zlib.h>
#include <iostream>
#include <fstream>


// [[Rcpp::export]]
int filterCells(
    std::string fragments,
    std::string outfile,
    Rcpp::Nullable<Rcpp::StringVector> keep_cells = R_NilValue,
    bool verbose = true
) {
  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(fragments.c_str(), "rb");

  // open new output file
  // TODO this creates an uncompressed file
  // ideally the output would also be bgzipped
  std::ofstream ofile;
  ofile.open(outfile);

  // return 1 if it can't find the file
  if (fileHandler == NULL) {
    Rcpp::Rcerr << "can't open file" << std::flush;
    return 1;
  }

  // C based buffered string parsing
  char* cb_char;
  size_t line_counter {1};
  uint32_t buffer_length = 256;
  char *buffer = new char[buffer_length];

  // Hash Map storing the barcodes to keep
  std::unordered_map<std::string, size_t> index_hash;

  size_t num_whitelist_cells {0};
  if (keep_cells.isNotNull()) {
    Rcpp::StringVector whitelist_cells(keep_cells);
    for (size_t i=0; i<whitelist_cells.size(); i++) {
      index_hash[Rcpp::as<std::string>(whitelist_cells[i])] = i;
    }
    num_whitelist_cells = index_hash.size();
    if (verbose) {
      Rcpp::Rcerr << "Keeping " << num_whitelist_cells
                  << " cell barcodes"
                  << std::endl << std::flush;
    }
  }

  // char * to string extraction
  std::string cb_seq;
  cb_seq.reserve(32);

  // looping over the fragments file
  while(gzgets(fileHandler, buffer, buffer_length) !=0 ){
      cb_char = strtok ( buffer, "\t" );

    for (auto i=1; i<=4; i++) {
      cb_char = strtok (NULL, "\t");
      switch(i) {
        case 3:
        // get cell barcode
        cb_seq.clear();
        cb_seq.append(cb_char);
        break;
        }
    }

    // if cell barcode in the keep list, write whole line to a new file
    auto cellcount = index_hash.count(cb_seq);
    if (cellcount > 0) {
      // cell present in requested set
      // write buffer to new file
      ofile << buffer;
      // TODO buffer is not what we want here
    }

    line_counter += 1;
    if (verbose) {
      if (line_counter % 10000000 == 0) {
        Rcpp::Rcerr << "\r                                                  ";
      }

      if (line_counter % 1000000 == 0) {
        Rcpp::Rcerr << "\rDone Processing " << line_counter / 1000000
                    << " million lines";
      }
    }
    if (line_counter % 10000000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  //Cleanup
  gzclose(fileHandler);
  ofile.close();

  return 0;
}
