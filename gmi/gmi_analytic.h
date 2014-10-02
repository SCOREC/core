/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_ANALYTIC_H
#define GMI_ANALYTIC_H

/** \file gmi_analytic.h
  \brief GMI analytic model interface */

#include "gmi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief the analytic parameterization of a model boundary
  \param p the input parametric coordinates. see gmi_eval for the
           format of this array
  \param x the resuting 3D point in space */
typedef void (*gmi_analytic_fun)(double const p[2], double x[3]);

/** \brief make an empty analytic model */
struct gmi_model* gmi_make_analytic(void);
/** \brief add an entity to the analytic model
  \param m the analytic model
  \param dim the dimension of the entity
  \param tag the dimension-unique tag for the entity
  \param f the analytic function defining the entity.
           set to NULL for interior entities.
  \param periodic for each dimension (d) of the entity,
                  gmi_periodic will return periodic[d]
  \param ranges for each dimension (d) of the entity,
                gmi_range will return ranges[d] */
void gmi_add_analytic(struct gmi_model* m, int dim, int tag,
    gmi_analytic_fun f, int* periodic, double (*ranges)[2]);
void gmi_add_analytic_region(struct gmi_model* m, int tag);

#ifdef __cplusplus
}
#endif

#endif

