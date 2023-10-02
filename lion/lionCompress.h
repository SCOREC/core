#ifndef LION_COMPRESS
#define LION_COMPRESS

namespace lion {

extern const bool can_compress;

void compress(void* dest, unsigned long& destLen,
    const void* source, unsigned long sourceLen);

unsigned long compressBound(unsigned long sourceLen);

}

#endif
