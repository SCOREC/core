/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maSize.h"
#include "apfMatrix.h"
#include <apfShape.h>
#include <cstdlib>
#include <pcu_util.h>

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

static void orthogonalizeR(Matrix& R)
{
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

static void orthogonalEigenDecompForSymmetricMatrix(Matrix const& A, Vector& v, Matrix& R)
{
  /* here we assume A to be real symmetric 3x3 matrix, 
   * we should be able to get 3 orthogonal eigen vectors
   * we also normalize the eigen vectors */
  double eigenValues[3];
  Vector eigenVectors[3];
  
  apf::eigen(A, eigenVectors, eigenValues);
  
  Matrix RT(eigenVectors); // eigen vectors are stored in the rows of RT

  RT[0] = RT[0].normalize();
  RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
  RT[1] = RT[1].normalize();
  RT[2] = apf::cross(RT[0],RT[1]);

  v = Vector(eigenValues);
  R = transpose(RT);
}

/* the length, area, or volume of
   the parent element for this
   entity type */
static double parentMeasure[apf::Mesh::TYPES] =
{0.0     //vert
,2.0     //edge
,1.0/2.0 //tri
,4.0     //quad
,1.0/6.0 //tet
,-42.0   //hex - definitely not sure
,-42.0   //prism - definitely not sure
,8.0/3.0 //pyramid
};

class SizeFieldIntegrator : public apf::Integrator
{
  public:
    SizeFieldIntegrator(SizeField* sF):
      Integrator(1),
      measurement(0),
      sizeField(sF),
      meshElement(0),
      dimension(0)
    {
    }
    void process(apf::MeshElement* me)
    {
      this->inElement(me);
      int np = countIntPoints(me,this->order);
      for (int p=0; p < np; ++p)
      {
        ipnode = p;
        Vector point;
        getIntPoint(me,this->order,p,point);
        double w = getIntWeight(me,this->order,p);
        this->atPoint(point,w,0);
      }
      this->outElement();
    }
    void inElement(apf::MeshElement* me)
    {
      meshElement = me;
      dimension = apf::getDimension(me);
    }
    void atPoint(Vector const& p , double w, double )
    {
      Matrix Q;
      sizeField->getTransform(meshElement,p,Q);
      Matrix J;
      apf::getJacobian(meshElement,p,J);
/* transforms the rows of J, the differential tangent vectors,
   into the metric space, then uses the generalized determinant */
      double dV2 = apf::getJacobianDeterminant(J*Q,dimension);
      measurement += w*dV2;
    }
    double measurement;
    SizeField* sizeField;
  private:
    apf::MeshElement* meshElement;
    int dimension;
};

struct MetricSizeField : public SizeField
{
  double measure(Entity* e)
  {
    SizeFieldIntegrator sFI(this); 
    apf::MeshElement* me = apf::createMeshElement(mesh, e);
    sFI.process(me);
    apf::destroyMeshElement(me);
    return sFI.measurement;
  }
  bool shouldSplit(Entity* edge)
  {
    return this->measure(edge) > 1.5;
  }
  bool shouldCollapse(Entity* edge)
  {
    return this->measure(edge) < 0.5;
  }
  double getWeight(Entity* e)
  {
    /* parentMeasure is used to normalize */
    return measure(e) / parentMeasure[mesh->getType(e)];
  }
  Mesh* mesh;
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
  BothEval()
  {
  }
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
  SizesEval()
  {
  }
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
  FrameEval()
  {
  }
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

struct LogMEval : public apf::Function
{
  LogMEval()
  {
  }
  LogMEval(AnisotropicFunction* f)
  {
    function = f;
    cachedVert = 0;
  }
  void updateCache(Entity* v)
  {
    if (v == cachedVert)
      return;
    Matrix R;
    Vector h;
    function->getValue(v, R, h);
    Matrix S( -2*log(h[0]),0,0,
              0,-2*log(h[1]),0,
              0,0,-2*log(h[2]));
    cachedLogM = R*S*transpose(R);
    cachedVert = v;
  }
  void getLogM(Entity* v, Matrix& f)
  {
    updateCache(v);
    f = cachedLogM;
  }
  void eval(Entity* e, double* result)
  {
    Matrix* f = (Matrix*) result;
    getLogM(e, *f);
  }
  Entity* cachedVert;
  Matrix cachedLogM;
  AnisotropicFunction* function;
};

struct AnisoSizeField : public MetricSizeField
{
  AnisoSizeField()
  {
  }
  AnisoSizeField(Mesh* m, AnisotropicFunction* f):
    bothEval(f),
    sizesEval(&bothEval),
    frameEval(&bothEval)
  {
    mesh = m;
    hField = apf::createUserField(m, "ma_sizes", apf::VECTOR,
        apf::getLagrange(1), &sizesEval);
    rField = apf::createUserField(m, "ma_frame", apf::MATRIX,
        apf::getLagrange(1), &frameEval);
  }
  ~AnisoSizeField()
  {
    apf::destroyField(hField);
    apf::destroyField(rField);
  }
  void init(Mesh* m, apf::Field* sizes, apf::Field* frames)
  {
    mesh = m;
    hField = sizes;
    rField = frames;
  }
  void getTransform(
      apf::MeshElement* me,
      Vector const& xi,
      Matrix& Q)
  {
    apf::Element* hElement = apf::createElement(hField,me);
    apf::Element* rElement = apf::createElement(rField,me);
    Vector h;
    Matrix R;
    apf::getVector(hElement,xi,h);
    apf::getMatrix(rElement,xi,R);
    apf::destroyElement(hElement);
    apf::destroyElement(rElement);
    orthogonalizeR(R);
    Matrix S(1/h[0],0,0,
             0,1/h[1],0,
             0,0,1/h[2]);
    Q = R*S;
  }
  void interpolate(
      apf::MeshElement* parent,
      Vector const& xi,
      Entity* newVert)
  {
    apf::Element* rElement = apf::createElement(rField,parent);
    apf::Element* hElement = apf::createElement(hField,parent);
    Vector h;
    apf::getVector(hElement,xi,h);
    Matrix R;
    apf::getMatrix(rElement,xi,R);
    orthogonalizeR(R);
    this->setValue(newVert,R,h);
    apf::destroyElement(hElement);
    apf::destroyElement(rElement);
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
  apf::Field* hField;
  apf::Field* rField;
  BothEval bothEval;
  SizesEval sizesEval;
  FrameEval frameEval;
};

struct LogAnisoSizeField : public MetricSizeField
{
  LogAnisoSizeField(Mesh* m, AnisotropicFunction* f):
    logMEval(f)
  {
    mesh = m;
    logMField = apf::createUserField(m, "ma_logM", apf::MATRIX,
        apf::getLagrange(1), &logMEval);
  }
  ~LogAnisoSizeField()
  {
    apf::destroyField(logMField);
  }
  void init(Mesh* m, apf::Field* logM)
  {
    mesh = m;
    logMField = logM;
  }
  void getTransform(
      apf::MeshElement* me,
      Vector const& xi,
      Matrix& Q)
  {
    apf::Element* logMElement = apf::createElement(logMField,me);
    Matrix logM;
    apf::getMatrix(logMElement,xi,logM);
    apf::destroyElement(logMElement);
    Vector v;
    Matrix R;
    orthogonalEigenDecompForSymmetricMatrix(logM, v, R);
    Matrix S( sqrt(exp(v[0])), 0, 0,
              0, sqrt(exp(v[1])), 0,
              0, 0, sqrt(exp(v[2])));
    Q = R*S;
  }
  void interpolate(
      apf::MeshElement* parent,
      Vector const& xi,
      Entity* newVert)
  {
    apf::Element* logMElement = apf::createElement(logMField,parent);
    Matrix logM;
    apf::getMatrix(logMElement,xi,logM);
    this->setValue(newVert,logM);
    apf::destroyElement(logMElement);
  }
  void setValue(
      Entity* vert,
      Matrix const& logM)
  {
    apf::setMatrix(logMField,vert,0,logM);
  }
  void setIsotropicValue(
      Entity* vert,
      double value)
  {
    this->setValue(vert,
                   Matrix(value,0,0,
                          0,value,0,
                          0,0,value));
  }
  apf::Field* logMField;
  LogMEval logMEval;
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
      PCU_ALWAYS_ASSERT(apf::getValueType(field)==apf::SCALAR);
      PCU_ALWAYS_ASSERT(apf::getShape(field)==apf::getLagrange(1));
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
  AnisoSizeField* f = new AnisoSizeField();
  f->init(m, sizes, frames);
  return f;
}

SizeField* makeSizeField(Mesh* m, AnisotropicFunction* f, int const interpolationOption)
{
  // interpolationOption is 0 by default
  if(interpolationOption == 0) 
    return new AnisoSizeField(m, f);
  else
    return new LogAnisoSizeField(m,f);
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
