/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_NULL_H
#define GMI_NULL_H

#include "gmi.h"

#ifdef __cplusplus
extern "C" {
#endif

void gmi_register_null(void);

struct gmi_ent* gmi_null_find(struct gmi_model* m, int dim, int tag);

#ifdef __cplusplus
}
#endif

#endif
