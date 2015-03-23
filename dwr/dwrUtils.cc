/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <iostream>
#include <PCU.h>
#include "dwrUtils.h"

namespace dwr {

void print(const char* format, ...)
{
  if (PCU_Comm_Self())
    return;
  printf("DWR: ");
  va_list ap;
  va_start(ap,format);
  vfprintf(stdout,format,ap);
  va_end(ap);
  printf("\n");
}

}
