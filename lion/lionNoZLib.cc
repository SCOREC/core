#include "lionCompress.h"

#include <cstdlib>

namespace lion {

const bool can_compress = false;

void compress(void* dest, unsigned long& destLen,
    const void* source, unsigned long sourceLen)
{
  (void) dest;
  (void) destLen;
  (void) source;
  (void) sourceLen;
  abort();
}

}

