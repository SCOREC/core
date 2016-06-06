/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <iostream>
#include <vector>
#include <mpi.h>
#include "apf.h"
#include "apfMDS.h"

#include "apfOmega_h.h"

pMesh pumi_ghost_create (int brgType, int ghostType, int numLayer, int includeCopy)
{
  
  if (brgType!=0 && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: bridge type "<<brgType<<" not supported\n";
    return NULL;
  }
  
  if (ghostType!=pumi::instance()->mesh->getDimension() && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" faild: ghost type "<<brgType<<" not supported\n";
    return NULL;
  }
  
  if (numLayer<1 && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid numLayer "<<numLayer<<"\n";
    return NULL;
  }

  if (includeCopy==0 && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: includeCopy=0"<<" not supported\n";
    return NULL;
  }
  // back up the original mesh for recovery
  if (!pumi::instance()->org_mesh) 
  {
    pumi::instance()->org_mesh = pumi::instance()->mesh;

    // build orig_node_flag for "isghost" query
    pumi::instance()->org_node_flag = new std::vector<bool>;
    int num_global_vtx;
    MPI_Allreduce(&(pumi::instance()->num_own_vtx), &num_global_vtx, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
    pumi::instance()->org_node_flag->resize(num_global_vtx);

    for (std::vector<bool>::iterator vit=pumi::instance()->org_node_flag->begin();
          vit!=pumi::instance()->org_node_flag->end(); ++vit)
      *vit = false;

    apf::MeshEntity* e;
    apf::MeshIterator* it = pumi::instance()->org_mesh->begin(0);
    while ((e = pumi::instance()->org_mesh->iterate(it)))
      pumi::instance()->org_node_flag->at(pumi_ment_getglobalid(e)) = true;
    pumi::instance()->org_mesh->end(it);
  }
  else // re-ghosting
  {
    if (!pumi_rank()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: accumulative ghosting not supported\n";
      return NULL;
  }

  // create ghosted mesh
  pMesh gm = apf::makeEmptyMdsMesh(pumi::instance()->model, ghostType, false);
  osh_t osh_mesh = osh::fromAPF(pumi::instance()->mesh);
  osh_ghost(osh_mesh, numLayer);
  osh::toAPF(osh_mesh, gm);
  osh_free(osh_mesh);

  pumi::instance()->mesh = gm;
  return gm;
}

void pumi_ghost_delete (pMesh m)
{
  if (!pumi_rank()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: not supported\n";
  return;
}

void pumi_ghost_info (pMesh m, std::vector<int>& ghostinfo)
{
  if (!pumi_rank()) 
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: not supported\n";
}

