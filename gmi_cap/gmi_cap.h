/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_CAP_H
#define GMI_CAP_H


#include "gmi.h"
#include "CapstoneModule.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Analysis.h"
#include "CreateMG_Framework_Application.h"
#include "CreateMG_Framework_Attributes.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Geometry.h"
#include "CreateMG_Framework_Mesh.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

gmi_ent* toGmiEntity(M_GTopo topo);
M_GTopo fromGmiEntity(gmi_ent* g);

void gmi_cap_start(void);
void gmi_cap_stop(void);
void gmi_register_cap(void);

/* struct gmi_model* gmi_cap_load(const char* nativefile, const char* smdfile); */
struct gmi_model* gmi_import_cap(GeometryDatabaseInterface* gi);
GeometryDatabaseInterface* gmi_export_cap(struct gmi_model* m);

#endif


