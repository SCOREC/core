/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef PUMI_H
#define PUMI_H
#include <gmi.h>
#include <apfMesh2.h>

typedef apf::Mesh2* pMesh;
typedef gmi_model* pGeom;

//************************************
//************************************
//      0- SYSTEM-LEVEL FUNCTIONS
//************************************
//************************************
void pumi_start();
void pumi_finalize(bool do_mpi_finalize=true);

int pumi_size();
int pumi_rank();

void pumi_sync(void);
void pumi_info();

//************************************
//************************************
//      1- MESH FUNCTIONS
//************************************
//************************************

//************************************
//  1.1 Mesh Instance Management 
//************************************
pGeom pumi_geom_create(const char* fileName, const char* model_type="mesh");
pMesh pumi_mesh_create(pGeom geom, const char* fileName, int option, const char* mesh_type="mds");
void pumi_mesh_write (pMesh mesh, const char* fileName, const char* mesh_type="mds");
void pumi_mesh_delete(pMesh mesh);
#endif
