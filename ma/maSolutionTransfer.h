/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SOLUTIONTRANSFER_H
#define MA_SOLUTIONTRANSFER_H

#include <apf.h>
#include "maMesh.h"

namespace ma {

class SolutionTransfer
{
  public:
    virtual ~SolutionTransfer();
    virtual bool hasNodesOn(int dimension) = 0;
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert);
    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities);
    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities);
    int getTransferDimension();
};

SolutionTransfer* createFieldTransfer(apf::Field* f);

class SolutionTransfers : public SolutionTransfer
{
  public:
    SolutionTransfers();
    virtual ~SolutionTransfers();
    void add(SolutionTransfer* t);
    virtual bool hasNodesOn(int dimension);
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert);
    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities);
    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities);
  private:
    typedef std::vector<SolutionTransfer*> Transfers;
    Transfers transfers;
};

class AutoSolutionTransfer : public SolutionTransfers
{
  public:
    AutoSolutionTransfer(Mesh* m);
};

}

#endif
