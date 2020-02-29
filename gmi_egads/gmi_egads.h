/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_EGADS_H
#define GMI_EGADS_H

/** \file gmi_egads.h
  \brief GMI EGADS model interface */

#ifdef __cplusplus
extern "C" {
#endif

struct egObject;

/** \brief start the EGADS session */
void gmi_egads_start(void);
/** \brief end the EGADS session */
void gmi_egads_stop(void);
/** \brief register the EGADS model reader for .egads files */
void gmi_register_egads(void);

/** \brief load an EGADS file into a gmi_model object */
struct gmi_model* gmi_egads_load(const char* filename);

/** \brief initialize a gmi_model with an EGADS body and number of regions */
struct gmi_model* gmi_egads_init(struct egObject *body, int numRegions);

/** \brief initialize the model adjacency table for 3D regions */
void gmi_egads_init_adjacency(int ***adjacency);

#ifdef __cplusplus
}
#endif

#endif
