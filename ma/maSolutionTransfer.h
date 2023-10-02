/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SOLUTIONTRANSFER_H
#define MA_SOLUTIONTRANSFER_H

/** \file maSolutionTransfer.h
  \brief MeshAdapt solution transfer interface */

#include <apf.h>
#include "maMesh.h"

namespace ma {

/** \brief user-defined solution transfer base
 \details one of these classes should be defined for
          each field the user wants transferred during
          adaptivity. Then, use ma::SolutionTransfers
          to put them together into a single object
          to give to the MeshAdapt driver functions. */
class SolutionTransfer
{
  public:
    /** \brief user-defined destructor */
    virtual ~SolutionTransfer();
    /** \brief return true if this field has
               nodes on entities of this dimension */
    virtual bool hasNodesOn(int dimension) = 0;
    /** \brief perform solution transfer on refined vertex
      \details this vertex was created while refining an
               element. Given coordinates in the parent
               element's parametric space, transfer the solution
               from the parent element to the vertex
      \param parent MeshElement for the parent.
                    can be given to apf::createElement
                    so that element interpolation can be used.
      \param xi Coordinates of the vertex in the parent
                element's parametric space */
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert);
    /** \brief perform solution transfer on refined entities
      \details when there are nodes on entities other
               than vertices, it becomes necessary to transfer
               solution to these as well.
               If this is not the case, don't override this.
      \param parent the parent element
      \param newEntities entities created during refinement.
                         this may contain some entities that
                         don't have nodes, so check them */
    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities);
    /** \brief perform solution transfer on cavity replacement
      \details this will be called for operations like edge
      collapses and edge swaps. Both the old and new cavity
      are given, and the code should transfer solution into
      into the new cavity.
      \param oldElements all elements in the old cavity
      \param newEntities newly created entities which need transfer
      */
    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities);
    /** \brief for internal MeshAdapt use */
    int getTransferDimension();
};

/** \brief Creates a default solution transfer object for a field
  \details MeshAdapt has good algorithms for transferring
  nodal fields as well as using the Voronoi system for transferring
  integration point fields. */
SolutionTransfer* createFieldTransfer(apf::Field* f);

/** \brief a meta-object that carries out a series of transfers
  \details use this class to put together solution transfer
  objects for several fields before giving them to MeshAdapt. */
class SolutionTransfers : public SolutionTransfer
{
  public:
    SolutionTransfers();
    virtual ~SolutionTransfers();
    /** \brief add another field transfer object */
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

/** \brief MeshAdapt's automatic solution transfer system.
  \details will call ma::createFieldTransfer on all fields associated
  with the mesh and put them together. */
class AutoSolutionTransfer : public SolutionTransfers
{
  public:
    AutoSolutionTransfer(Mesh* m);
};

}

#endif
