#include <PCU.h>
#include "apfMesh2.h"
#include "apfField.h"
#include "apf.h"

#include "apfShape.h"
#include "apfTagData.h"

namespace apf
{

void Mesh2::setPoint(MeshEntity* e, int node, Vector3 const& p)
{
  setVector(Mesh::coordinateField,e,node,p);
}

MeshEntity* Mesh2::createVertex(ModelEntity* c, Vector3 const& point,
    Vector3 const& param)
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
    Vector3 const* points)
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

static void packTagClone(Mesh2* m, MeshTag* t, int to)
{
  std::string name;
  name = m->getTagName(t);
  packString(name, to);
  int type;
  type = m->getTagType(t);
  PCU_COMM_PACK(to, type);
  int size;
  size = m->getTagSize(t);
  PCU_COMM_PACK(to, size);
}

static MeshTag* unpackTagClone(Mesh2* m)
{
  std::string name = unpackString();
  int type;
  PCU_COMM_UNPACK(type);
  int size;
  PCU_COMM_UNPACK(size);
  if (type == apf::Mesh::DOUBLE)
    return m->createDoubleTag(name.c_str(), size);
  if (type == apf::Mesh::INT)
    return m->createIntTag(name.c_str(), size);
  if (type == apf::Mesh::LONG)
    return m->createLongTag(name.c_str(), size);
  return 0;
}

static void packTagClones(Mesh2* m, int to)
{
  DynamicArray<MeshTag*> tags;
  m->getTags(tags);
  int n = tags.getSize();
  PCU_COMM_PACK(to, n);
  /* warning! this loop goes backward to cater to MDS
     implementation-specific behavior.
     please forgive me. */
  for (int i = n - 1; i >= 0; --i)
    packTagClone(m, tags[i], to);
}

static void unpackTagClones(Mesh2* m)
{
  int n;
  PCU_COMM_UNPACK(n);
  for (int i = 0; i < n; ++i)
    unpackTagClone(m);
}

static void packFieldClone(Field* f, int to)
{
  std::string name = f->getName();
  packString(name, to);
  int valueType = f->getValueType();
  PCU_COMM_PACK(to, valueType);
  int components = f->countComponents();
  PCU_COMM_PACK(to, components);
  std::string shapeName = f->getShape()->getName();
  packString(shapeName, to);
  /* warning! this only supports tag-stored fields */
}

static Field* unpackFieldClone(Mesh2* m)
{
  std::string name = unpackString();
  int valueType;
  PCU_COMM_UNPACK(valueType);
  int components;
  PCU_COMM_UNPACK(components);
  std::string shapeName = unpackString();
  FieldShape* shape = getShapeByName(shapeName.c_str());
  /* warning! this only supports tag-stored fields */
  return makeField(m, name.c_str(), valueType, components, shape, new TagDataOf<double>);
}

static void packFieldClones(Mesh2* m, int to)
{
  int n = m->countFields();
  PCU_COMM_PACK(to, n);
  for (int i = 0; i < n; ++i)
    packFieldClone(m->getField(i), to);
}

static void unpackFieldClones(Mesh2* m)
{
  int n;
  PCU_COMM_UNPACK(n);
  for (int i = 0; i < n; ++i)
    unpackFieldClone(m);
}

void packDataClone(Mesh2* m, int to)
{
  packTagClones(m, to);
  packFieldClones(m, to);
}

void unpackDataClone(Mesh2* m)
{
  unpackTagClones(m);
  unpackFieldClones(m);
}

}
