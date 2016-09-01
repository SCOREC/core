/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"

pGeomTag pumi_geom_createTag (pGeom model, const char* tagName, int tagType, int tagSize);
void pumi_geom_deleteTag (pGeom model, pGeomTag tag, int forceDel);
pGeomTag pumi_geom_findTag (pGeom model, const char* tagName);
bool pumi_geom_hasTag (pGeom model, const pGeomTag tag);
int pumi_geom_getTag (pGeom model, std::vector<pGeomTag>& tags);
bool pumi_geom_isTagInUse (pGeom model, pGeomTag tag);

bool pumi_gent_hasTag (pGeomEnt ent, pGeomTag tag);
void pumi_gent_deleteTag (pGeomEnt ent, pGeomTag tag);
void pumi_gent_getTag (pGeomEnt ent, std::vector<pGeomTag>& tags);

//void pumi_gent_setPtrTag (pGeomEnt ent, pGeomTag tag, void* data);
//int pumi_gent_getPtrTag (pGeomEnt ent, pGeomTag tag, void** data);
void pumi_gent_setIntTag (pGeomEnt ent, pGeomTag tag, const int data);
int pumi_gent_getIntTag (pGeomEnt ent, pGeomTag tag);
void pumi_gent_setLongTag (pGeomEnt ent, pGeomTag tag, const long data);
long pumi_gent_getLongTag (pGeomEnt ent, pGeomTag tag);
void pumi_gent_setDblTag (pGeomEnt ent, pGeomTag tag, const double data);
double pumi_gent_getDblTag (pGeomEnt ent, pGeomTag tag);
//int pumi_gent_setEntTag (pGeomEnt ent, pGeomTag tag, const pGeomEnt data);
//int pumi_gent_getEntTag (pGeomEnt ent, pGeomTag tag, pGeomEnt* data);

void pumi_gent_setIntArrTag (pGeomEnt ent, pGeomTag tag, const int* data);
void pumi_gent_getIntArrTag (pGeomEnt ent, pGeomTag tag, int** data, int* data_size);
void pumi_gent_setLongArrTag (pGeomEnt ent, pGeomTag tag, const long* data);
void pumi_gent_getLongArrTag (pGeomEnt ent, pGeomTag tag, long** data, int* data_size);
void pumi_gent_setDblArrTag (pGeomEnt ent, pGeomTag tag, const double* data);
void pumi_gent_getDblArrTag (pGeomEnt ent, pGeomTag tag, double** data, int* data_size);
//void pumi_gent_setEntArrTag (pGeomEnt ent, pGeomTag tag, const pGeomEnt* data);
//int pumi_gent_getEntArrTag (pGeomEnt ent, pGeomTag tag, pGeomEnt** data, int* data_size);
//void pumi_gent_setStringTag(pGeomEnt ent, pGeomTag tag, const char* s);
//int pumi_gent_getStringTag(pGeomEnt ent, pGeomTag tag, const char*& s);
