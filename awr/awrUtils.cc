/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrUtils.h"
#include <PCU.h>
#include <stdlib.h>

namespace awr {

void fail(const char* why)
{
  if (!PCU_Comm_Self())
    fprintf(stderr,"AWR FAILED: %s\n",why);
  abort();
}

}
