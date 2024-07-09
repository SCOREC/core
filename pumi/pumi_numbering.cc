/****************************************************************************** 

  (c) 2004-2018 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include "apf.h"
#include "apfShape.h"
#include "apfFieldData.h"
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include <assert.h>
#include <iostream>

//************************************
// Node numbering
//************************************

pNumbering pumi_numbering_create
   (pMesh m, const char* name, pShape shape, int num_component)
{
  pNumbering n = m->findNumbering(name);
  if (n) 
  {
    if (!pumi_rank()) 
      std::cout<<"[PUMI INFO] "<<__func__<<" failed: numbering \""<<name<<"\" already exists\n";
    return n;
  }

  if (!shape) shape= m->getShape();
  return createNumbering(m, name, shape, num_component);
}

pNumbering pumi_numbering_createLocal (pMesh m, const char* name, pShape shape)
{
  pNumbering n = m->findNumbering(name);
  if (n) 
  {
    if (!pumi_rank()) 
      std::cout<<"[PUMI INFO] "<<__func__<<" failed: numbering \""<<name<<"\" already exists\n";
    return n;
  }

  if (!shape) shape= m->getShape();
  return numberOverlapNodes(m, name, shape);
}

pNumbering pumi_numbering_createGlobal(pMesh m, const char* name, pShape s, pOwnership o)
{
  pNumbering n = m->findNumbering(name);
  if (n) 
  {
    if (!pumi_rank()) 
      std::cout<<"[PUMI INFO] "<<__func__<<" failed: numbering \""<<name<<"\" already exists\n";
    return n;
  }

  if (!s) s= m->getShape();
  n = numberOwnedNodes(m, name, s, o);
  apf::globalize(n, m->getPCU());
  apf::synchronizeFieldData<int>(n->getData(), o, false); //synchronize(n, o); 
  return n;
}

pNumbering pumi_numbering_createOwn (pMesh m, const char* name, pShape shape, pOwnership o)
{
  pNumbering n = m->findNumbering(name);
  if (n) 
  {
    if (!pumi_rank()) 
      std::cout<<"[PUMI INFO] "<<__func__<<" failed: numbering \""<<name<<"\" already exists\n";
    return n;
  }

   if (!shape) shape= m->getShape();
   return numberOwnedNodes(m, name, shape, o);
}

pNumbering pumi_numbering_createOwnDim (pMesh m, const char* name, int dim, pOwnership o)
{  
  pNumbering n = m->findNumbering(name);
  if (n) 
  {
    if (!pumi_rank()) 
      std::cout<<"[PUMI INFO] "<<__func__<<" failed: numbering \""<<name<<"\" already exists\n";
    return n;
  }

  return numberOwnedDimension(m, name, dim, o);
}

// requirement: a mesh entity on process group "p" is owned by a process within group "p"
pNumbering pumi_numbering_createProcGrp (
    pMesh m, const char* name, int num_proc_grp, 
    int dim, pOwnership o)
{
  assert(m->getPCU()->Peers()%num_proc_grp==0);

  pNumbering n = m->findNumbering(name);
  if (n) 
  {
    if (!pumi_rank()) 
      std::cout<<"[PUMI INFO] "<<__func__<<" failed: numbering \""<<name<<"\" already exists\n";
    return n;
  }

  int self = m->getPCU()->Self();
  int pgrp_size = m->getPCU()->Peers()/num_proc_grp;
  int local_pgrpid = self/pgrp_size; // divide
  int pgrp_rank = self % pgrp_size;

  bool delete_ownership=false;
  if (!o) 
  {
    o = getSharing(m);
    delete_ownership=true;
  }

  //n = numberOwnedDimension(m, name, dim, o);  
  n = apf::createNumbering(m,name,apf::getConstant(dim),1);

  int owned_node_cnt=0;
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
  {
    if (!o->isOwned(e))
      continue;
    owned_node_cnt += n->countNodesOn(e);
  }

  int* in = new int;
  int* out_arr = new int[m->getPCU()->Peers()]; // out[i] has local_numOwnedPartBdryEnt of process i on all processes
  *in = owned_node_cnt;

  MPI_Allgather(in, 1, MPI_INT, out_arr, 1, MPI_INT, m->getPCU()->GetMPIComm());

  it = m->begin(dim);
  int nbr = 0;
  for (int pid=1; pid<=pgrp_rank; ++pid)
    nbr += out_arr[local_pgrpid*pgrp_size+pid-1];

  while ((e = m->iterate(it)))
  {
    if (!o->isOwned(e))
      continue;
    int nnodes = n->countNodesOn(e);
    for (int node=0; node < nnodes; ++node)
       number(n,e,node,0, nbr++);
  }
  m->end(it);

  apf::synchronizeFieldData<int>(n->getData(), o, false);
  if (delete_ownership) 
    delete o;
  delete [] out_arr;
  return n;
}


void pumi_numbering_delete(pNumbering n)
{
  destroyNumbering(n);
}

int pumi_numbering_getNumNode(pNumbering n)
{
  return apf::countNodes(n);
}

void pumi_node_setNumber(pNumbering nb, pMeshEnt e, int n, int c, int number)
{
  apf::number(nb, e, n, c, number);
}

int pumi_node_getNumber(pNumbering nb, pMeshEnt e, int n, int c)
{
  return apf::getNumber(nb, e, n, c);
}

bool pumi_node_isNumbered(pNumbering nb, pMeshEnt e, int n, int c)
{
  return apf::isNumbered(nb, e, n, c);
}

void pumi_numbering_print(pNumbering n, int pid)
{
  pMesh m = pumi::instance()->mesh;
  pMeshEnt e;
  int nnodes;
  pShape s = n->getShape();
  if (pid==-1) pid = m->getPCU()->Peers();
  for (int rank=0; rank<pid; ++rank)
  {
    if (rank==m->getPCU()->Self())
    {
      for(int dd = 0; dd < m->getDimension(); ++dd)
      {
        if (s->hasNodesIn(dd))
        {
          pMeshIter it = m->begin(dd);
          while((e = m->iterate(it)))
          {
            nnodes = n->countNodesOn(e);
            for (int node=0; node < nnodes; ++node)
               std::cout<<"("<<m->getPCU()->Self()<<") ent "<<pumi_ment_getID(e)
                     <<"/node "<<node<<" # "<<getNumber(n,e,node,0)<<" (ghost? "<<(m->isGhost(e)?1:0)<<")\n";
          } // while
          m->end(it);
        } // if shp
      } // for dd
    }  // if (pid==PCU_Comm_Self())
    MPI_Barrier(MPI_COMM_WORLD);
  } // for
}
