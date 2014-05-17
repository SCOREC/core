/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SIZE_H
#define MA_SIZE_H

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

class IdentitySizeField : public SizeField
{
  public:
    IdentitySizeField(Mesh* m);
    virtual double measure(Entity* e);
    virtual bool shouldSplit(Entity* edge);
    virtual bool shouldCollapse(Entity* edge);
    virtual double placeSplit(Entity* edge);
    virtual void interpolate(
        apf::MeshElement* parent,
        Vector const& xi,
        Entity* newVert);
    virtual void getTransform(
        apf::MeshElement* e,
        Vector const& xi,
        Matrix& t);
    virtual double getWeight(Entity* e);
    Mesh* getMesh() {return mesh;}
  private:
    Mesh* mesh;
};

class UniformRefiner : public IdentitySizeField
{
  public:
    UniformRefiner(Mesh* m);
    virtual bool shouldSplit(Entity* edge);
};

class MetricSizeField : public IdentitySizeField
{
  public:
    MetricSizeField(
        Mesh* m,
        std::string const& name = "ma_size");
    ~MetricSizeField();
    virtual double measure(Entity* e);
    virtual bool shouldSplit(Entity* edge);
    virtual bool shouldCollapse(Entity* edge);
    virtual double placeSplit(Entity* edge);
    virtual void interpolate(
        apf::MeshElement* parent,
        Vector const& xi,
        Entity* newVert);
    virtual void getTransform(
        apf::MeshElement* e,
        Vector const& xi,
        Matrix& t);
    virtual double getWeight(Entity* e);
    void setValue(
        Entity* vert,
        Matrix const& r,
        Vector const& h);
    void setIsotropicValue(
        Entity* vert,
        double value);
    apf::Field* rField;
    apf::Field* hField;
};

class AnisotropicFunction
{
  public:
    virtual ~AnisotropicFunction();
    virtual void getValue(Entity* vert, Matrix& r, Vector& h) = 0;
};

void initialize(MetricSizeField* field, AnisotropicFunction* function);

class IsotropicFunction
{
  public:
    virtual ~IsotropicFunction();
    virtual double getValue(Entity* vert) = 0;
};

void initialize(MetricSizeField* field, IsotropicFunction* function);

double getAverageEdgeLength(Mesh* m);

void printIsotropic(MetricSizeField* sf);
void printAnisotropic(MetricSizeField* sf);

}

#endif
