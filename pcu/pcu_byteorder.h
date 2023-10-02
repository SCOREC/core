// derived from
// http://stackoverflow.com/questions/2100331/c-macro-definition-to-determine-big-endian-or-little-endian-machine
#ifndef PCU_BYTEORDER_H
#define PCU_BYTEORDER_H

#include <limits.h>
#include <stdint.h>

#if CHAR_BIT != 8
#error "unsupported char size"
#endif

enum {
  PCU_LITTLE_ENDIAN = 0x03020100ul,
  PCU_BIG_ENDIAN = 0x00010203ul
};

static const union {
  unsigned char bytes[4];
  uint32_t value;
} pcu_host_order = { { 0, 1, 2, 3 } };

#define PCU_HOST_ORDER (pcu_host_order.value)

#endif
