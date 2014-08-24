/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_NULL_H
#define GMI_NULL_H

/** \file gmi_null.h
  \brief GMI null model interface */

#include "gmi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief register null model system
  \details
  The null model system is the shame of GMI.
  It exists to create small meshes by hand for debugging and
  to deal with meshes created by less than intelligent methods.
  Please avoid it at all costs.

  The null model system registers with file
  extension ".null", so any filename that ends with that
  will create an empty null model without reading any real files.

  gmi_find calls to the null model will return an existing pointer
  if found, otherwise a new entity will be created with that
  tag and dimension and that new pointer will be returned.
  This allows a classified mesh to build up a set of model entities
  while being loaded, or users can just manually call
  model find functions to "create" new entities. */
void gmi_register_null(void);

/** \brief the null model implementation of gmi_find
  \details this is not for public users, just call gmi_find */
struct gmi_ent* gmi_null_find(struct gmi_model* m, int dim, int tag);

#ifdef __cplusplus
}
#endif

#endif
