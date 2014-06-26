/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_ANALYTIC_H
#define GMI_ANALYTIC_H

#include "gmi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*gmi_analytic_fun)(double const p[2], double x[3]);

struct gmi_model* gmi_make_analytic(void);
void gmi_add_analytic(struct gmi_model* m, int dim, int tag,
    gmi_analytic_fun f, int* periodic, double (*ranges)[2]);

#ifdef __cplusplus
}
#endif

#endif

