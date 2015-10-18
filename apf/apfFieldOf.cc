#include "apfFieldOf.h"
#include "apfElementOf.h"

namespace apf {

template <class T>
class Project : public FieldOp
{
  public:
    Project(FieldOf<T>* a, FieldOf<T>* b)
    {
      to = a;
      from = b;
      mesh = a->getMesh();
    }
    bool inEntity(MeshEntity* e)
    {
      meshElement = createMeshElement(mesh,e);
      fromElement = static_cast<ElementOf<T>*>(createElement(from,meshElement));
      return true;
    }
    void atNode(int n)
    {
      Vector3 xi;
      to->getShape()->getNodeXi(fromElement->getType(),n,xi);
      T value[1];
      value[0] = fromElement->getValue(xi);
      to->setNodeValue(fromElement->getEntity(),n,value);
    }
    void outEntity()
    {
      destroyMeshElement(meshElement);
      destroyElement(fromElement);
    }
    void run()
    {
      apply(to);
    }
    FieldOf<T>* to;
    FieldOf<T>* from;
    Mesh* mesh;
    VectorElement* meshElement;
    ElementOf<T>* fromElement;
};

template <class T>
void project(FieldOf<T>* to, FieldOf<T>* from)
{
  Project<T> projector(to,from);
  projector.run();
}

/* instantiate here */
template void project<double>(FieldOf<double>* to, FieldOf<double>* from);
template void project<Vector3>(FieldOf<Vector3>* to, FieldOf<Vector3>* from);
template void project<Matrix3x3>(FieldOf<Matrix3x3>* to, FieldOf<Matrix3x3>* from);

template <class T>
class Axpy : public FieldOp
{
  public:
    void run(double a_, FieldOf<T>* x_, FieldOf<T>* y_)
    {
      a = a_;
      x = x_;
      y = y_;
      apply(y);
    }
    bool inEntity(MeshEntity* e)
    {
      entity = e;
      return true;
    }
    void atNode(int n)
    {
      T xv[1];
      T yv[1];
      x->getNodeValue(entity, n, xv);
      y->getNodeValue(entity, n, yv);
      T axpyv[1];
      axpyv[0] = (xv[0] * a) + yv[0];
      y->setNodeValue(entity, n, axpyv);
    }
    double a;
    FieldOf<T>* x;
    FieldOf<T>* y;
    MeshEntity* entity;
};

template <class T>
void axpy(double a, FieldOf<T>* x, FieldOf<T>* y)
{
  Axpy<T> op;
  op.run(a,x,y);
}

/* instantiate here */
template void axpy<double>(double, FieldOf<double>*, FieldOf<double>*);
template void axpy<Vector3>(double, FieldOf<Vector3>*, FieldOf<Vector3>*);
template void axpy<Matrix3x3>(double, FieldOf<Matrix3x3>*, FieldOf<Matrix3x3>*);

}
