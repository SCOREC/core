#include <pcu_util.h>
#include "apfMesh2.h"
#include "apfField.h"
#include "apf.h"
#include "apfShape.h"
#include "apfTagData.h"
#include "apfNumbering.h"

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
    BuildCallback* cb,
    bool* p_made)
{
  MeshEntity* e = findUpward(m,type,down);
  if (e) {
    if (p_made) *p_made = false;
    return e;
  }
  e = m->createEntity(type,c,down);
  if (cb) cb->call(e);
  if (p_made) *p_made = true;
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
  m->getPCU()->Pack(to,t);
  m->getPCU()->Pack(to,e);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  m->getPCU()->Pack(to,nd);
  for (int i = 0; i < nd; ++i) {
    Copies remotes;
    m->getRemotes(down[i], remotes);
    MeshEntity* dr = remotes[to];
    m->getPCU()->Pack(to,dr);
  }
}

static void unpackProposal(int& t, MeshEntity*& e, Downward& da, pcu::PCU *PCUObj)
{
  PCUObj->Unpack(t);
  PCUObj->Unpack(e);
  int nd;
  PCUObj->Unpack(nd);
  for (int i=0; i < nd; ++i)
    PCUObj->Unpack(da[i]);
}

void stitchMesh(Mesh2* m)
{
  initResidence(m, 0);
  int d_max = m->getDimension();
  MeshEntity* e;
  for (int d=1; d < d_max; ++d) {
    m->getPCU()->Begin();
    MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it))) {
      Parts candidateParts;
      getCandidateParts(m, e, candidateParts);
      candidateParts.erase(m->getId());
      APF_ITERATE(Parts, candidateParts, pit)
        packProposal(m, e,*pit);
    }
    m->end(it);
    m->getPCU()->Send();
    while (m->getPCU()->Listen()) {
      int from = m->getPCU()->Sender();
      while (!m->getPCU()->Unpacked()) {
        int t;
        Downward da;
        unpackProposal(t, e, da, m->getPCU());
        MeshEntity* found = findUpward(m, t, da);
        if (found)
          m->addRemote(found, from, e);
      }
    }
    initResidence(m, d);
  }
}

static void packTagClone(Mesh2* m, MeshTag* t, int to, pcu::PCU *PCUObj)
{
  std::string name;
  name = m->getTagName(t);
  packString(name, to, PCUObj);
  int type;
  type = m->getTagType(t);
  PCUObj->Pack(to, type);
  int size;
  size = m->getTagSize(t);
  PCUObj->Pack(to, size);
}

static MeshTag* unpackTagClone(Mesh2* m)
{
  std::string name = unpackString(m->getPCU());
  int type;
  m->getPCU()->Unpack(type);
  int size;
  m->getPCU()->Unpack(size);
  if (type == apf::Mesh::DOUBLE)
    return m->createDoubleTag(name.c_str(), size);
  if (type == apf::Mesh::INT)
    return m->createIntTag(name.c_str(), size);
  if (type == apf::Mesh::LONG)
    return m->createLongTag(name.c_str(), size);
  return 0;
}

static void packTagClones(Mesh2* m, int to, pcu::PCU *PCUObj)
{
  DynamicArray<MeshTag*> tags;
  m->getTags(tags);
  int n = tags.getSize();
  PCUObj->Pack(to, n);
  
  /* warning! this loop goes backward to cater to MDS
     implementation-specific behavior.
     please forgive me. */
  for (int i = n - 1; i >= 0; --i)
    packTagClone(m, tags[i], to, PCUObj);
}

static void unpackTagClones(Mesh2* m)
{
  int n;
  m->getPCU()->Unpack(n);
  for (int i = 0; i < n; ++i)
    unpackTagClone(m);
}

static void packFieldClone(Field* f, int to, pcu::PCU *PCUObj)
{
  std::string name = f->getName();
  packString(name, to, PCUObj);
  int valueType = f->getValueType();
  PCUObj->Pack(to, valueType);
  
  int components = f->countComponents();
  PCUObj->Pack(to, components);
  
  std::string shapeName = f->getShape()->getName();
  packString(shapeName, to, PCUObj);

  /* warning! this only supports tag-stored fields */
}

static Field* unpackFieldClone(Mesh2* m)
{
  std::string name = unpackString(m->getPCU());
  int valueType;
  m->getPCU()->Unpack(valueType);
  int components;
  m->getPCU()->Unpack(components);
  std::string shapeName = unpackString(m->getPCU());
  FieldShape* shape = getShapeByName(shapeName.c_str());
  PCU_ALWAYS_ASSERT(shape);
  /* warning! this only supports tag-stored fields */
  return makeField(m, name.c_str(), valueType, components, shape, new TagDataOf<double>);
}

static void packFieldClones(Mesh2* m, int to, pcu::PCU *PCUObj)
{
  int n = m->countFields();
  PCUObj->Pack(to, n);
  
  for (int i = 0; i < n; ++i)
    packFieldClone(m->getField(i), to, PCUObj);
}

static void unpackFieldClones(Mesh2* m)
{
  int n;
  m->getPCU()->Unpack(n);
  for (int i = 0; i < n; ++i)
    unpackFieldClone(m);
}

static void packMeshShape(Mesh2* m, int to, pcu::PCU *PCUObj)
{
  std::string shapeName = m->getShape()->getName();
  packString(shapeName, to, PCUObj);
}

static void unpackMeshShape(Mesh2* m)
{
  std::string shapeName = unpackString(m->getPCU());
  FieldShape* shape = getShapeByName(shapeName.c_str());
  PCU_ALWAYS_ASSERT(shape);
  if (shape != m->getShape()) {
    m->changeShape(shape, /*project=*/false);
  }
}

void packDataClone(Mesh2* m, int to, pcu::PCU *PCUObj)
{
  packTagClones(m, to, PCUObj);
  packMeshShape(m, to, PCUObj);
  packFieldClones(m, to, PCUObj);
}

void unpackDataClone(Mesh2* m)
{
  unpackTagClones(m);
  unpackMeshShape(m);
  unpackFieldClones(m);
}

void clear(Mesh2* m) {
  while (m->countFields()) destroyField(m->getField(0));
  while (m->countNumberings()) destroyNumbering(m->getNumbering(0));
  while (m->countGlobalNumberings()) {
    destroyGlobalNumbering(m->getGlobalNumbering(0));
  }
  DynamicArray<MeshTag*> tags;
  m->getTags(tags);
  for (size_t i = 0; i < tags.getSize(); ++i) {
    m->destroyTag(tags[i]);
  }
  m->clear_();
}

}
