/******************************************************************************

  Copyright 2011 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "pcu_util.h"
#include <stdlib.h>

void PCU_Assert_Fail(const char* msg) {
  fprintf(stderr, "%s", msg);
  abort();
}
