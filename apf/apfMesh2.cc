#include "apfMesh2.h"
#include "apfField.h"
#include "apf.h"

namespace apf
{

void Mesh2::setPoint(MeshEntity* e, int node, Vector3 const& p)
{
  setVector(Mesh::coordinateField,e,node,p);
}

MeshEntity* Mesh2::createVertex(ModelEntity* c, Vector3 const& point, Vector3 const& param)
{
  MeshEntity* v = createVert(c);
  setPoint(v,0,point);
  setParam(v,param);
  return v;
}

void displaceMesh(Mesh2* m, Field* d, double factor)
{
  m->getCoordinateField()->axpy(factor,d);
}

MeshEntity* makeOrFind(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** down,
    BuildCallback* cb)
{
  MeshEntity* e = findUpward(m,type,down);
  if (e) return e;
  e = m->createEntity(type,c,down);
  if (cb) cb->call(e);
  return e;
}

class ElementBuilder : public ElementVertOp
{
  public:
    ElementBuilder(
        Mesh2* m,
        ModelEntity* c,
        BuildCallback* cb)
    {
      mesh = m;
      modelEntity = c;
      callback = cb;
    }
    virtual MeshEntity* apply(int type, MeshEntity** down)
    {
      return makeOrFind(mesh,modelEntity,type,down,callback);
    }
  private:
    Mesh2* mesh;
    ModelEntity* modelEntity;
    BuildCallback* callback;
};

MeshEntity* buildElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    MeshEntity** verts,
    BuildCallback* cb)
{
  ElementBuilder b(m,c,cb);
  return b.run(type,verts);
}

MeshEntity* buildOneElement(
    Mesh2* m,
    ModelEntity* c,
    int type,
    Vector3* points)
{
  int nv = Mesh::adjacentCount[type][0];
  Downward v;
  for (int i=0; i < nv; ++i)
    v[i] = m->createVertex(c,points[i],Vector3(0,0,0));
  return buildElement(m,c,type,v);
}

}
