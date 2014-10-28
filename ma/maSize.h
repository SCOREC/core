/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SIZE_H
#define MA_SIZE_H

/** \file maSize.h
    \brief MeshAdapt Size Fields */

#include <apf.h>
#include "maMesh.h"

namespace ma {

typedef apf::Matrix3x3 Matrix;

class SizeField
{
  public:
    virtual ~SizeField();
    virtual double measure(Entity* e) = 0;
    virtual bool shouldSplit(Entity* edge) = 0;
    virtual bool shouldCollapse(Entity* edge) = 0;
    virtual double placeSplit(Entity* edge) = 0;
    virtual void interpolate(
        apf::MeshElement* parent,
        Vector const& xi,
        Entity* newVert) = 0;
    virtual void getTransform(
        apf::MeshElement* e,
        Vector const& xi,
        Matrix& t) = 0;
    virtual double getWeight(Entity* e) = 0;
};

struct IdentitySizeField : public SizeField
{
  IdentitySizeField(Mesh* m);
  double measure(Entity* e);
  bool shouldSplit(Entity*);
  bool shouldCollapse(Entity*);
  double placeSplit(Entity*);
  void interpolate(
      apf::MeshElement* parent,
      Vector const& xi,
      Entity* newVert);
  void getTransform(
          apf::MeshElement*,
          Vector const&,
          Matrix& t);
  double getWeight(Entity*);
  Mesh* mesh;
};

struct UniformRefiner : public IdentitySizeField
{
  UniformRefiner(Mesh* m):
    IdentitySizeField(m)
  {
  }
  bool shouldSplit(Entity*)
  {
    return true;
  }
};

/** \brief User-defined Anisotropic size function */
class AnisotropicFunction
{
  public:
    virtual ~AnisotropicFunction();
    /** \brief get the size field value at this vertex
      \param r the orthonormal basis frame
      \param h the desired element sizes along each
               of the frame's basis vectors */
    virtual void getValue(Entity* vert, Matrix& r, Vector& h) = 0;
};

/** \brief User-defined Isotropic size function */
class IsotropicFunction
{
  public:
    virtual ~IsotropicFunction();
    /** \brief get the desired element size at this vertex */
    virtual double getValue(Entity* vert) = 0;
};

SizeField* makeSizeField(Mesh* m, apf::Field* sizes, apf::Field* frames);
SizeField* makeSizeField(Mesh* m, apf::Field* size);
SizeField* makeSizeField(Mesh* m, AnisotropicFunction* f);
SizeField* makeSizeField(Mesh* m, IsotropicFunction* f);

double getAverageEdgeLength(Mesh* m);

}

#endif
