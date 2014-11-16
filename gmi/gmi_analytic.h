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

#include "gmi_base.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief the analytic parameterization of a model boundary
  \param p the input parametric coordinates. see gmi_eval for the
           format of this array
  \param x the resuting 3D point in space
  \param u pointer to user data */
typedef void (*gmi_analytic_fun)(double const p[2], double x[3], void* u);

/** \brief a re-parametrization from one entity to another
  \param from the coordinates in the input parametric space
  \param to   the coordinates in the output parametric space
  \param u    extra user-provided data */
typedef void (*gmi_reparam_fun)(double const from[2], double to[2], void* u);

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
                gmi_range will return ranges[d]
  \param user_data pointer that will be passed to analytic
                   function for this entity */
struct gmi_ent* gmi_add_analytic(struct gmi_model* m, int dim, int tag,
    gmi_analytic_fun f, int* periodic, double (*ranges)[2], void* user_data);
/** \brief get the analytic user data
  \details this function returns the pointer passed as (user_data)
  to gmi_add_analytic when creating entity (e) */
void* gmi_analytic_data(struct gmi_model* m, struct gmi_ent* e);

#define gmi_analytic_topo gmi_base_topo

/** \brief add a re-parameterization to the model
  \details given a use (u) which is a piece of a boundary,
  define the map between the parametric space of the used
  entity and the parametric space of the using entity.
  \param m the analytic model
  \param u the model use (see agm.h for how to build this)
  \param f the reparameterization code
  \param user_data extra data for the re-parameterization code
  */
void gmi_add_analytic_reparam(struct gmi_model* m, struct agm_use u,
    gmi_reparam_fun f, void* user_data);

/** \brief get the re-parameterization user data
  \details this function returns the pointer passed as (user_data)
  to gmi_add_analytic_reparam when creating entity (e) */
void* gmi_analytic_reparam_data(struct gmi_model* m, struct agm_use u);

/** \brief create a non-parametric model entity
  \details this creates a highest-dimensional model entity.
           typically users do not create a parametric definition
           of the 3D regions or 2D faces, since this is equivalent
           to the analytic function. */
void gmi_add_analytic_cell(struct gmi_model* m, int dim, int tag);
/** \brief gmi_add_analytic_cell(m, 3, tag)
  \details this is provided as a bit of backwards compatibility */
void gmi_add_analytic_region(struct gmi_model* m, int tag);

#ifdef __cplusplus
}
#endif

#endif

