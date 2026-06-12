#ifndef SIGNAC_FRAGMENT_UTILS_H
#define SIGNAC_FRAGMENT_UTILS_H

#include <zlib.h>
#include <cstring>

// Check whether the buffer just filled by gzgets() holds a complete line.
// gzgets() reads at most buffer_length - 1 characters, stopping early only at a
// newline or end-of-file. A complete line therefore either ends in a newline,
// or coincides with the end of the file (the final line may have no trailing
// newline). If neither holds, the line was longer than the buffer and has been
// split mid-record; parsing it would silently corrupt the data.
inline bool fragmentLineComplete(const char* buffer, gzFile file) {
  size_t n = std::strlen(buffer);
  if (n > 0 && buffer[n - 1] == '\n') {
    return true;
  }
  return gzeof(file) != 0;
}

#endif  // SIGNAC_FRAGMENT_UTILS_H
