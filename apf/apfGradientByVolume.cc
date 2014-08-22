#include <apfCavityOp.h>
#include <apf.h>

namespace apf {

template<class T>
class GradientOf;

template<>
class GradientOf<double>
{
  public:
    typedef Vector3 type;
};

template<>
class GradientOf<Vector3>
{
  public:
    typedef Matrix3x3 type;
};

template<class T>
class GetGradient;

template<>
class GetGradient<double>
{
  public:
    Vector3 operator()(Element* e, Vector3 const& xi)
    {
      Vector3 v;
      getGrad(e,xi,v);
      return v;
    }
};

template<>
class GetGradient<Vector3>
{
  public:
    Matrix3x3 operator()(Element* e, Vector3 const& xi)
    {
      Matrix3x3 m;
      getVectorGrad(e,xi,m);
      return m;
    }
};

template <class T>
class GradientIntegrator : public Integrator
{
  public:
    typedef typename GradientOf<T>::type GT;
    GradientIntegrator(Field* f_in):Integrator(1)
              /* assuming all linear !  -------^ */
    {
      f = f_in;
      isFirst = true;
      volumeSum = 0;
    }
    virtual void inElement(MeshElement* me)
    {
      e = createElement(f,me);
    }
    virtual void atPoint(Vector3 const& p, double w, double dV)
    {
      GT grad = getGradient(e,p);
      if (isFirst)
      {
        gradSum = grad*w*dV;
        isFirst = false;
      }
      else
        gradSum = gradSum + grad*w*dV;
      volumeSum += w*dV;
    }
    virtual void outElement()
    {
      destroyElement(e);
    }
    GT getResult() {return gradSum / volumeSum;}
  private:
    Field* f;
    Element* e;
    bool isFirst;
    GT gradSum;
    double volumeSum;
    GetGradient<T> getGradient;
};

template<class T>
class SetValue;

template<>
class SetValue<Vector3>
{
  public:
    void operator()(Field* f, MeshEntity* v, Vector3 const& value)
    {
      setVector(f,v,0,value);
    }
};

template<>
class SetValue<Matrix3x3>
{
  public:
    void operator()(Field* f, MeshEntity* v, Matrix3x3 const& value)
    {
      setMatrix(f,v,0,value);
    }
};

template<class T>
class RecoverGradient : public CavityOp
{
  public:
    typedef typename GradientOf<T>::type GT;
    RecoverGradient(Field* f_in, Field* gradf_in):
      CavityOp(getMesh(f_in))
    {
      mesh = getMesh(f_in);
      f = f_in;
      gradf = gradf_in;
    }
    virtual Outcome setEntity(MeshEntity* v)
    {
      if (hasEntity(gradf,v))
        return SKIP;
      if ( ! requestLocality(&v,1))
        return REQUEST;
      vert = v;
      return OK;
    }
    virtual void apply()
    {
      Adjacent elements;
      mesh->getAdjacent(vert,mesh->getDimension(),elements);
      GradientIntegrator<T> integrator(f);
      for (size_t i=0; i < elements.getSize(); ++i)
      {
        MeshEntity* e = elements[i];
        MeshElement* me = createMeshElement(mesh,e);
        integrator.process(me);
        destroyMeshElement(me);
      }
      GT grad = integrator.getResult();
      setValue(gradf,vert,grad);
    }
  private:
    Mesh* mesh;
    MeshEntity* vert;
    Field* f;
    Field* gradf;
    SetValue<GT> setValue;
};

template<class T>
void recoverGradientByVolume(Field* f, Field* gradf)
{
  RecoverGradient<T> op(f,gradf);
  op.applyToDimension(0);
}

Field* recoverGradientByVolume(Field* f)
{
  Mesh* mesh = getMesh(f);
  std::string name = getName(f);
  name += "_vol";
  int valueType = getValueType(f);
  Field* gradf = createLagrangeField(mesh,name.c_str(),valueType+1,1);
  if (valueType==SCALAR)
    recoverGradientByVolume<double>(f,gradf);
  else
  { assert(valueType==VECTOR);
    recoverGradientByVolume<Vector3>(f,gradf);
  }
  return gradf;
}

}
