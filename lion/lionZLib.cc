#include "lionCompress.h"

#include <zlib.h>

namespace lion {

const bool can_compress = true;

void compress(void* dest, unsigned long& destLen,
    const void* source, unsigned long sourceLen)
{
  ::compress((Bytef*)dest, &destLen, (const Bytef*)source, sourceLen);
}

}
