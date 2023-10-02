/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef GMI_LOOKUP_H
#define GMI_LOOKUP_H

#include "agm.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmi_lookup;

struct gmi_lookup* gmi_new_lookup(struct agm* topo);
void gmi_free_lookup(struct gmi_lookup* l);
void gmi_set_lookup(struct gmi_lookup* l, struct agm_ent e, int tag);
int gmi_get_lookup(struct gmi_lookup* l, struct agm_ent e);
struct agm_ent gmi_look_up(struct gmi_lookup* l, enum agm_ent_type t, int tag);
void gmi_freeze_lookup(struct gmi_lookup* l, enum agm_ent_type t);
void gmi_unfreeze_lookups(struct gmi_lookup* l);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
