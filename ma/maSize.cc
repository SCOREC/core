/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSize.h"
#include <cstdlib>
#include <PCU.h>

namespace ma {

SizeField::~SizeField()
{
}

IdentitySizeField::IdentitySizeField(Mesh* m):
  mesh(m)
{
}

double IdentitySizeField::measure(Entity* e)
{
  apf::MeshElement* me = apf::createMeshElement(mesh,e);
  double x = apf::measure(me);
  apf::destroyMeshElement(me);
  return x;
}

bool IdentitySizeField::shouldSplit(Entity*)
{
  return false;
}

bool IdentitySizeField::shouldCollapse(Entity*)
{
  return false;
}

double IdentitySizeField::placeSplit(Entity*)
{
  return 0.5;
}

void IdentitySizeField::interpolate(
        apf::MeshElement*,
        Vector const&,
        Entity*)
{
}

void IdentitySizeField::getTransform(
        apf::MeshElement*,
        Vector const&,
        Matrix& t)
{
  t = Matrix(1,0,0,
             0,1,0,
             0,0,1);
}

double IdentitySizeField::getWeight(Entity*)
{
  return 1.0;
}

UniformRefiner::UniformRefiner(Mesh* m):
  IdentitySizeField(m)
{
}

bool UniformRefiner::shouldSplit(Entity*)
{
  return true;
}

MetricSizeField::MetricSizeField(
    Mesh* m,
    std::string const& name):
  IdentitySizeField(m)
{
  std::string fieldName = name;
  fieldName += "_r";
  rField = apf::createLagrangeField(m,fieldName.c_str(),apf::MATRIX,1);
  fieldName = name;
  fieldName += "_h";
  hField = apf::createLagrangeField(m,fieldName.c_str(),apf::VECTOR,1);
}

MetricSizeField::~MetricSizeField()
{
  apf::destroyField(rField);
  apf::destroyField(hField);
}

void makeQ(Matrix const& R, Vector const& h, Matrix& Q)
{
  Matrix S(1/h[0],0,0,
           0,1/h[1],0,
           0,0,1/h[2]);
  Q = R*S;
}

void interpolateR(
    apf::Element* rElement,
    Vector const& xi,
    Matrix& R)
{
  apf::getMatrix(rElement,xi,R);
  Matrix RT = transpose(R);
  /* we store R and h separately because the interpolation
     is not simple: we re-normalize the rotation matrix R
     after interpolating its components. */
  for (int i=0; i < 3; ++i)
    RT[i] = RT[i].normalize();
  R = transpose(RT);
}

void interpolateQ(
    apf::Element* rElement,
    apf::Element* hElement,
    Vector const& xi,
    Matrix& Q)
{
  Matrix R;
  interpolateR(rElement,xi,R);
  Vector h;
  apf::getVector(hElement,xi,h);
  makeQ(R,h,Q);
}

class SizeFieldIntegrator : public apf::Integrator
{
  public:
    SizeFieldIntegrator(MetricSizeField* f):
      Integrator(1),
      measurement(0),
      sizeField(f)
    {}
    virtual void inElement(apf::MeshElement* me)
    {
      dimension = apf::getDimension(me);
      rElement = apf::createElement(sizeField->rField,me);
      hElement = apf::createElement(sizeField->hField,me);
    }
    virtual void outElement()
    {
      apf::destroyElement(rElement);
      apf::destroyElement(hElement);
    }
    virtual void atPoint(Vector const& p, double w, double)
    {
      Matrix Q;
      interpolateQ(rElement,hElement,p,Q);
      Matrix J;
      apf::getJacobian(apf::getMeshElement(rElement),p,J);
/* transforms the rows of J, the differential tangent vectors,
   into the metric space, then uses the generalized determinant */
      double dV2 = apf::getJacobianDeterminant(J*Q,dimension);
      measurement += w*dV2;
    }
    double measurement;
  private:
    MetricSizeField* sizeField;
    int dimension;
    apf::Element* rElement;
    apf::Element* hElement;
};

double MetricSizeField::measure(Entity* e)
{
  SizeFieldIntegrator integrator(this);
  apf::MeshElement* me = apf::createMeshElement(getMesh(),e);
  integrator.process(me);
  apf::destroyMeshElement(me);
  return integrator.measurement;
}

bool MetricSizeField::shouldSplit(Entity* edge)
{
  return this->measure(edge) > 1.5;
}

bool MetricSizeField::shouldCollapse(Entity* edge)
{
  return this->measure(edge) < 0.5;
}

double MetricSizeField::placeSplit(Entity* edge)
{
  Mesh* m = getMesh();
  Entity* v[2];
  m->getDownward(edge,0,v);
  Vector p[2];
  m->getPoint(v[0],0,p[0]);
  m->getPoint(v[1],0,p[1]);
  Vector e = (p[1]-p[0]).normalize();
  double h[2];
  for (int i=0; i < 2; ++i)
  {
    Matrix R;
    apf::getMatrix(rField,v[i],0,R);
    Vector hv;
    apf::getVector(hField,v[i],0,hv);
    Matrix Q;
    makeQ(R,hv,Q);
    Vector e2 = transpose(Q)*e;
    h[i] = 1/(e2.getLength());
  }
/* the approximation given by Li based on the assumption
   that desired length varies linearly along the edge. */
  return 1.0/(1.0+sqrt(h[1]/h[0]));
}

void MetricSizeField::interpolate(
    apf::MeshElement* parent,
    Vector const& xi,
    Entity* newVert)
{
  apf::Element* rElement = apf::createElement(rField,parent);
  apf::Element* hElement = apf::createElement(hField,parent);
  Matrix R;
  interpolateR(rElement,xi,R);
  Vector h;
  apf::getVector(hElement,xi,h);
  this->setValue(newVert,R,h);
  apf::destroyElement(hElement);
  apf::destroyElement(rElement);
}

void MetricSizeField::getTransform(
        apf::MeshElement* e,
        Vector const& xi,
        Matrix& t)
{
  apf::Element* rElement = apf::createElement(rField,e);
  apf::Element* hElement = apf::createElement(hField,e);
  interpolateQ(rElement,hElement,xi,t);
  apf::destroyElement(rElement);
  apf::destroyElement(hElement);
}

/* the length, area, or volume of
   the parent element for this
   entity type */
double parentMeasure[TYPES] =
{0.0     //vert
,2.0     //edge - not sure
,1.0/2.0 //tri
,4.0     //quad - not sure
,1.0/6.0 //tet
,-42.0   //hex - definitely not sure
,-42.0   //prism - definitely not sure
,-42.0   //pyramid - definitely not sure
};

double MetricSizeField::getWeight(Entity* e)
{
  Mesh* m = getMesh();
  /* parentMeasure is used to normalize */
  return measure(e) / parentMeasure[m->getType(e)];
}

void MetricSizeField::setValue(
    Entity* vert,
    Matrix const& r,
    Vector const& h)
{
  apf::setMatrix(rField,vert,0,r);
  apf::setVector(hField,vert,0,h);
}

void MetricSizeField::setIsotropicValue(
    Entity* vert,
    double value)
{
  this->setValue(vert,
                 Matrix(1,0,0,
                        0,1,0,
                        0,0,1),
                 Vector(value,value,value));
}

AnisotropicFunction::~AnisotropicFunction()
{
}

void initialize(MetricSizeField* field, AnisotropicFunction* function)
{
  Mesh* m = field->getMesh();
  apf::MeshIterator* it = m->begin(0);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    Matrix r;
    Vector h;
    function->getValue(e,r,h);
    field->setValue(e,r,h);
  }
  m->end(it);
}

IsotropicFunction::~IsotropicFunction()
{
}

void initialize(MetricSizeField* field, IsotropicFunction* function)
{
  Mesh* m = field->getMesh();
  apf::MeshIterator* it = m->begin(0);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    double value = function->getValue(e);
    field->setIsotropicValue(e,value);
  }
  m->end(it);
}

double getAverageEdgeLength(Mesh* m)
{
  IdentitySizeField sizeField(m);
  double sums[2];
  double& length_sum = sums[0];
  double& edge_count = sums[1];
  length_sum = edge_count = 0;
  apf::MeshIterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    length_sum += sizeField.measure(e);
    edge_count += 1.0;
  }
  m->end(it);
  PCU_Add_Doubles(sums,2);
  return length_sum / edge_count;
}

void printIsotropic(MetricSizeField* sf)
{
  apf::Mesh* m = apf::getMesh(sf->rField);
  Tag* tag = m->createDoubleTag("size",1);
  Iterator* it = m->begin(0);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    Vector h;
    apf::getVector(sf->hField,e,0,h);
    double size = h[0];
    m->setDoubleTag(e,tag,&size);
  }
  m->end(0);
  apf::writeVtkFiles("size",m);
  apf::removeTagFromDimension(m,tag,0);
  m->destroyTag(tag);
}

void printAnisotropic(MetricSizeField* sf)
{
  apf::Mesh* m = apf::getMesh(sf->rField);
  Tag* tag_x = m->createDoubleTag("size_x",3);
  Tag* tag_y = m->createDoubleTag("size_y",3);
  Tag* tag_z = m->createDoubleTag("size_z",3);
  Iterator* it = m->begin(0);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    Vector h;
    apf::getVector(sf->hField,e,0,h);
    Matrix R;
    apf::getMatrix(sf->rField,e,0,R);
    Matrix Q;
    makeQ(R,h,Q);
    Matrix frame = apf::invert(Q);
    m->setDoubleTag(e,tag_x,&(frame[0][0]));
    m->setDoubleTag(e,tag_y,&(frame[1][0]));
    m->setDoubleTag(e,tag_z,&(frame[2][0]));
  }
  m->end(0);
  apf::writeVtkFiles("size",m);
  apf::removeTagFromDimension(m,tag_x,0);
  apf::removeTagFromDimension(m,tag_y,0);
  apf::removeTagFromDimension(m,tag_z,0);
  m->destroyTag(tag_x);
  m->destroyTag(tag_y);
  m->destroyTag(tag_z);
}

}
