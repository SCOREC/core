/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
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
#include <PCU.h>
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
  apf::globalize(n);
  apf::synchronizeFieldData<int>(n->getData(), o, false); //synchronize(n, o); 
  return n;
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


