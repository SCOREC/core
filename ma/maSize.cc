/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maSize.h"
#include <apfShape.h>
#include <cstdlib>

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

static void makeQ(Matrix const& R, Vector const& h, Matrix& Q)
{
  Matrix S(1/h[0],0,0,
           0,1/h[1],0,
           0,0,1/h[2]);
  Q = R*S;
}

static void interpolateR(
    apf::Element* rElement,
    Vector const& xi,
    Matrix& R)
{
  apf::getMatrix(rElement,xi,R);
  /* by the way, the principal direction vectors
     are in the columns, so lets do our cleanup
     work on the transpose */
  Matrix RT = transpose(R);
  /* we have to make sure this remains a rotation matrix.
     lets assume that component-wise interpolation gave
     a somewhat decent answer.
     (if not, the truly proper solution is to interpolate
      quaternion components, see also Lie groups)
     R[0] should be the interpolated main principal direction
     vector. lets first normalize it. */
  RT[0] = RT[0].normalize();
  /* now the next principal direction should roughly align
     with R[1], but first it must be orthogonal to R[0].
     so, remove its component along R[0] (remember R[0] is unit length)*/
  RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
  /* and normalize it too */
  RT[1] = RT[1].normalize();
  /* finally, forget what was in R[2] and just make it the cross
     product of R[0] and R[1] to make a nice frame */
  RT[2] = apf::cross(RT[0],RT[1]);
  /* and back out of the transpose */
  /* if you're wondering, this is why R and h are stored separately:
     so we can do this interpolation dance on the rotation matrix */
  R = transpose(RT);
}

static void interpolateQ(
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
    SizeFieldIntegrator(apf::Field* r, apf::Field* h):
      Integrator(1),
      measurement(0),
      rField(r),
      hField(h)
    {}
    virtual void inElement(apf::MeshElement* me)
    {
      dimension = apf::getDimension(me);
      rElement = apf::createElement(rField,me);
      hElement = apf::createElement(hField,me);
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
    apf::Field* rField;
    apf::Field* hField;
    int dimension;
    apf::Element* rElement;
    apf::Element* hElement;
};

/* the length, area, or volume of
   the parent element for this
   entity type */
static double parentMeasure[TYPES] =
{0.0     //vert
,2.0     //edge
,1.0/2.0 //tri
,4.0     //quad
,1.0/6.0 //tet
,-42.0   //hex - definitely not sure
,-42.0   //prism - definitely not sure
,8.0/3.0 //pyramid
};

struct MetricSizeField : public SizeField
{
  void init(Mesh* m, apf::Field* sizes, apf::Field* frames)
  {
    mesh = m;
    rField = frames;
    hField = sizes;
  }
  ~MetricSizeField()
  {
  }
  double measure(Entity* e)
  {
    SizeFieldIntegrator integrator(rField, hField);
    apf::MeshElement* me = apf::createMeshElement(mesh, e);
    integrator.process(me);
    apf::destroyMeshElement(me);
    return integrator.measurement;
  }
  bool shouldSplit(Entity* edge)
  {
    return this->measure(edge) > 1.5;
  }
  bool shouldCollapse(Entity* edge)
  {
    return this->measure(edge) < 0.5;
  }
  double placeSplit(Entity* edge)
  {
    Entity* v[2];
    mesh->getDownward(edge,0,v);
    Vector p[2];
    mesh->getPoint(v[0],0,p[0]);
    mesh->getPoint(v[1],0,p[1]);
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
  void interpolate(
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
  void getTransform(
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
  double getWeight(Entity* e)
  {
    /* parentMeasure is used to normalize */
    return measure(e) / parentMeasure[mesh->getType(e)];
  }
  void setValue(
      Entity* vert,
      Matrix const& r,
      Vector const& h)
  {
    apf::setMatrix(rField,vert,0,r);
    apf::setVector(hField,vert,0,h);
  }
  void setIsotropicValue(
      Entity* vert,
      double value)
  {
    this->setValue(vert,
                   Matrix(1,0,0,
                          0,1,0,
                          0,0,1),
                   Vector(value,value,value));
  }
  Mesh* mesh;
  apf::Field* rField;
  apf::Field* hField;
};

AnisotropicFunction::~AnisotropicFunction()
{
}

IsotropicFunction::~IsotropicFunction()
{
}

struct IsoWrapper : public AnisotropicFunction
{
  IsoWrapper(IsotropicFunction* f)
  {
    function = f;
  }
  void getValue(Entity* vert, Matrix& r, Vector& h)
  {
    r = Matrix(1,0,0,
               0,1,0,
               0,0,1);
    double s = function->getValue(vert);
    h = Vector(s,s,s);
  }
  IsotropicFunction* function;
};

struct BothEval
{
  BothEval(AnisotropicFunction* f)
  {
    function = f;
    cachedVert = 0;
  }
  void updateCache(Entity* v)
  {
    if (v == cachedVert)
      return;
    function->getValue(v, cachedFrame, cachedSizes);
    cachedVert = v;
  }
  void getSizes(Entity* v, Vector& s)
  {
    updateCache(v);
    s = cachedSizes;
  }
  void getFrame(Entity* v, Matrix& f)
  {
    updateCache(v);
    f = cachedFrame;
  }
  Entity* cachedVert;
  Vector cachedSizes;
  Matrix cachedFrame;
  AnisotropicFunction* function;
};

struct SizesEval : public apf::Function
{
  SizesEval(BothEval* b)
  {
    both = b;
  }
  void eval(Entity* e, double* result)
  {
    Vector* s = (Vector*) result;
    both->getSizes(e, *s);
  }
  BothEval* both;
};

struct FrameEval : public apf::Function
{
  FrameEval(BothEval* b)
  {
    both = b;
  }
  void eval(Entity* e, double* result)
  {
    Matrix* f = (Matrix*) result;
    both->getFrame(e, *f);
  }
  BothEval* both;
};

struct AnisoSizeField : public MetricSizeField
{
  AnisoSizeField(Mesh* m, AnisotropicFunction* f):
    bothEval(f),
    sizesEval(&bothEval),
    frameEval(&bothEval)
  {
    sizesField = apf::createUserField(m, "ma_sizes", apf::VECTOR,
        apf::getLagrange(1), &sizesEval);
    frameField = apf::createUserField(m, "ma_frame", apf::MATRIX,
        apf::getLagrange(1), &frameEval);
    this->init(m, sizesField, frameField);
  }
  ~AnisoSizeField()
  {
    apf::destroyField(sizesField);
    apf::destroyField(frameField);
  }
  BothEval bothEval;
  SizesEval sizesEval;
  FrameEval frameEval;
  apf::Field* sizesField;
  apf::Field* frameField;
};

struct IsoSizeField : public AnisoSizeField
{
  IsoSizeField(Mesh* m, IsotropicFunction* f):
    AnisoSizeField(m, &wrapper),
    wrapper(f)
  {
  }
  IsoWrapper wrapper;
};

class FieldReader : public IsotropicFunction
{
  public:
    FieldReader(apf::Field* f)
    {
      field = f;
      assert(apf::getValueType(field)==apf::SCALAR);
      assert(apf::getShape(field)==apf::getMesh(field)->getShape());
    }
    virtual ~FieldReader() {}
    virtual double getValue(Entity* vert)
    {
      return apf::getScalar(field,vert,0);
    }
    apf::Field* field;
};

struct IsoUserField : public IsoSizeField
{
  IsoUserField(Mesh* m, apf::Field* f):
    IsoSizeField(m, &reader),
    reader(f)
  {
  }
  FieldReader reader;
};

SizeField* makeSizeField(Mesh* m, apf::Field* sizes, apf::Field* frames)
{
  MetricSizeField* f = new MetricSizeField();
  f->init(m, sizes, frames);
  return f;
}

SizeField* makeSizeField(Mesh* m, AnisotropicFunction* f)
{
  return new AnisoSizeField(m, f);
}

SizeField* makeSizeField(Mesh* m, IsotropicFunction* f)
{
  return new IsoSizeField(m, f);
}

SizeField* makeSizeField(Mesh* m, apf::Field* size)
{
  return new IsoUserField(m, size);
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

}
