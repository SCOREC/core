#include "apfMesh2.h"
#include "apfField.h"
#include "apf.h"
#include <PCU.h>

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

void initResidence(apf::Mesh2* m, int d)
{
  apf::MeshIterator* it = m->begin(d);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::Copies remotes;
    m->getRemotes(e, remotes);
    apf::Parts parts;
    APF_ITERATE(apf::Copies, remotes, rit)
      parts.insert(rit->first);
    parts.insert(m->getId());
    m->setResidence(e, parts);
  }
  m->end(it);
}

static void intersect(
    Parts& a,
    Parts const& b)
{
  for (Parts::iterator it = a.begin(); it != a.end();) {
    if ( ! b.count(*it))
      a.erase(*(it++));
    else
      ++it;
  }
}

static void getCandidateParts(Mesh* m, MeshEntity* e, Parts& parts)
{
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  //down adjacencies have to at least be shared
  for (int i = 0; i < nd; ++i)
    if (!m->isShared(down[i]))
      return;
  bool first = true;
  for (int i = 0; i < nd; ++i) {
    MeshEntity* da = down[i];
    Parts rp;
    m->getResidence(da, rp);
    if (first) {
      parts = rp;
      first = false;
    } else {
      intersect(parts,rp);
    }
  }
}

static void packProposal(Mesh* m, MeshEntity* e, int to)
{
  int t = m->getType(e);
  PCU_COMM_PACK(to,t);
  PCU_COMM_PACK(to,e);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  PCU_COMM_PACK(to,nd);
  for (int i = 0; i < nd; ++i) {
    Copies remotes;
    m->getRemotes(down[i], remotes);
    MeshEntity* dr = remotes[to];
    PCU_COMM_PACK(to,dr);
  }
}

static void unpackProposal(int& t, MeshEntity*& e, Downward& da)
{
  PCU_COMM_UNPACK(t);
  PCU_COMM_UNPACK(e);
  int nd;
  PCU_COMM_UNPACK(nd);
  for (int i=0; i < nd; ++i)
    PCU_COMM_UNPACK(da[i]);
}

void stitchMesh(Mesh2* m)
{
  initResidence(m, 0);
  int d_max = m->getDimension();
  MeshEntity* e;
  for (int d=1; d < d_max; ++d) {
    PCU_Comm_Begin();
    MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it))) {
      Parts candidateParts;
      getCandidateParts(m, e, candidateParts);
      candidateParts.erase(m->getId());
      APF_ITERATE(Parts, candidateParts, pit)
        packProposal(m, e,*pit);
    }
    m->end(it);
    PCU_Comm_Send();
    while (PCU_Comm_Listen()) {
      int from = PCU_Comm_Sender();
      while (!PCU_Comm_Unpacked()) {
        int t;
        Downward da;
        unpackProposal(t, e, da);
        MeshEntity* found = findUpward(m, t, da);
        if (found)
          m->addRemote(found, from, e);
      }
    }
    initResidence(m, d);
  }
}

}
