/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "apfCoordData.h"
#include "apfVectorField.h"
#include "apfShape.h"
#include "apfNumbering.h"
#include "apfTagData.h"
#include <gmi.h>

namespace apf {

/* note: these tables must be consistent
   with those of the mesh database, and
   in the best design should be read from
   the database.
   The goal of supporting multiple databases
   makes this difficult for now */

int const tri_edge_verts[3][2] =
{{0,1},{1,2},{2,0}};

int const quad_edge_verts[4][2] =
{{0,1},{1,2},{2,3},{3,0}};

/*
   Canonical vertex and edge numbering
   for a tetrahedron, the middle (3) vertex
   is higher out of the page than the 0-1-2 triangle.
             2
            /|\
           / | \
          /  |  \
         /   |   \
        /    |5   \
      2/     |     \1
      /     _|_     \
     /   __/ 3 \__   \
    / _3/        4\__ \
   /_/               \_\
  0---------------------1
             0
*/

int const tet_edge_verts[6][2] =
{{0,1}
,{1,2}
,{2,0}
,{0,3}
,{1,3}
,{2,3}};

int const prism_edge_verts[9][2] =
{{0,1},{1,2},{2,0}
,{0,3},{1,4},{2,5}
,{3,4},{4,5},{5,3}};

int const pyramid_edge_verts[8][2] =
{{0,1},{1,2},{2,3},{3,0}
,{0,4},{1,4},{2,4},{3,4}};

int const tet_tri_verts[4][3] =
{{0,1,2}
,{0,1,3}
,{1,2,3}
,{0,2,3}};

int const prism_tri_verts[2][3] =
{{0,1,2}
,{3,4,5}};

int const prism_quad_verts[3][4] =
{{0,1,4,3}
,{1,2,5,4}
,{2,0,3,5}};

int const pyramid_tri_verts[4][3] =
{{0,1,4}
,{1,2,4}
,{2,3,4}
,{3,0,4}};

/* the following tables should be the same as
 * those found in the mesh databases, with the possible
 * exception of the behavior at the same dimension
 */

int const Mesh::adjacentCount[TYPES][4] = 
{{1,-1,-1,-1} //vertex
,{2,1,-1,-1} //edge
,{3,3,1,-1}//tri
,{4,4,1,-1}//quad
,{4,6,4,1}//tet
,{8,12,6,1}//hex
,{6,9,5,1}//prism
,{5,8,5,1}//pyramid
};

int const Mesh::typeDimension[TYPES] = 
{ 0, //vertex
  1, //edge
  2, //tri
  2, //quad
  3, //tet
  3, //hex
  3, //prism
  3 //pyramid
};

char const* const Mesh::typeName[TYPES] =
{"vertex",
 "edge",
 "triangle",
 "quad",
 "tet",
 "hex",
 "prism",
 "pyramid"
};

void Mesh::init(FieldShape* s)
{
  coordinateField = new VectorField();
  FieldBase* baseP = coordinateField;
  CoordData* data = new CoordData();
  baseP->init("coordinates",this,s,data);
  data->init(baseP);
  hasFrozenFields = false;
}

Mesh::~Mesh()
{
  delete coordinateField;
}

int Mesh::getModelType(ModelEntity* e)
{
  return gmi_dim(getModel(), (gmi_ent*)e);
}

int Mesh::getModelTag(ModelEntity* e)
{
  return gmi_tag(getModel(), (gmi_ent*)e);
}

ModelEntity* Mesh::findModelEntity(int type, int tag)
{
  return (ModelEntity*)gmi_find(getModel(), type, tag);
}

bool Mesh::canSnap()
{
  return gmi_can_eval(getModel());
}

void Mesh::snapToModel(ModelEntity* m, Vector3 const& p, Vector3& x)
{
  gmi_eval(getModel(), (gmi_ent*)m, &p[0], &x[0]);
}

void Mesh::getParamOn(ModelEntity* g, MeshEntity* e, Vector3& p)
{
  ModelEntity* from_g = toModel(e);
  if (g == from_g)
    return getParam(e, p);
  gmi_ent* from = (gmi_ent*)from_g;
  gmi_ent* to = (gmi_ent*)g;
  Vector3 from_p;
  getParam(e, from_p);
  gmi_reparam(getModel(), from, &from_p[0], to, &p[0]);
}

bool Mesh::getPeriodicRange(ModelEntity* g, int axis, double range[2])
{
  gmi_ent* e = (gmi_ent*)g;
  gmi_range(getModel(), e, axis, range);
  return gmi_periodic(getModel(), e, axis);
}

void Mesh::getPoint(MeshEntity* e, int node, Vector3& p)
{
  getVector(coordinateField,e,node,p);
}

FieldShape* Mesh::getShape() const
{
  return coordinateField->getShape();
}

void Mesh::changeCoordinateField(Field* f)
{
  delete coordinateField;
  coordinateField = f;
}

void Mesh::addField(Field* f)
{
  fields.push_back(f);
}

void Mesh::removeField(Field* f)
{
  fields.erase(std::find(fields.begin(),fields.end(),f));
}

Field* Mesh::findField(const char* name)
{
  std::string n(name);
  for (size_t i=0; i < fields.size(); ++i)
    if (n==getName(fields[i]))
      return fields[i];
  return 0;
}

int Mesh::countFields()
{
  return static_cast<int>(fields.size());
}

Field* Mesh::getField(int i)
{
  return fields[i];
}

void Mesh::addNumbering(Numbering* n)
{
  numberings.push_back(n);
}

void Mesh::removeNumbering(Numbering* n)
{
  numberings.erase(std::find(numberings.begin(),numberings.end(),n));
}

Numbering* Mesh::findNumbering(const char* name)
{
  std::string n(name);
  for (size_t i=0; i < numberings.size(); ++i)
    if (n==getName(numberings[i]))
      return numberings[i];
  return 0;
}

int Mesh::countNumberings()
{
  return static_cast<int>(numberings.size());
}

Numbering* Mesh::getNumbering(int i)
{
  return numberings[i];
}

void Mesh::addGlobalNumbering(GlobalNumbering* n)
{
  globalNumberings.push_back(n);
}

void Mesh::removeGlobalNumbering(GlobalNumbering* n)
{
  globalNumberings.erase(std::find(
        globalNumberings.begin(), globalNumberings.end(), n));
}

int Mesh::countGlobalNumberings()
{
  return static_cast<int>(globalNumberings.size());
}

GlobalNumbering* Mesh::getGlobalNumbering(int i)
{
  return globalNumberings[i];
}


void unite(Parts& into, Parts const& from)
{
  into.insert(from.begin(),from.end());
}

void getPeers(Mesh* m, int d, Parts& peers)
{
  assert(d < m->getDimension());
  MeshEntity* e;
  MeshIterator* ents = m->begin(d);
  while ((e = m->iterate(ents)))
  {
    Parts residence;
    m->getResidence(e,residence);
    unite(peers,residence);
  }
  m->end(ents);
  peers.erase(m->getId());
}

static bool residesOn(Mesh* m, MeshEntity* e, int part)
{
  Parts residence;
  m->getResidence(e,residence);
  return residence.count(part);
}

MeshEntity* iterateBoundary(Mesh* m, MeshIterator* it, int part)
{
  MeshEntity* e;
  while ((e = m->iterate(it))&&( ! residesOn(m,e,part)));
  return e;
}

Migration::Migration(Mesh* m)
{
  mesh = m;
  tag = m->createIntTag("apf_migrate",1);
}

Migration::Migration(Mesh* m, MeshTag* existingTag)
{
  mesh = m;
  tag = existingTag;
}

Migration::~Migration()
{
  for (size_t i=0; i < elements.size(); ++i)
    mesh->removeTag(elements[i],tag);
  mesh->destroyTag(tag);
}

int Migration::count()
{
  return elements.size();
}

MeshEntity* Migration::get(int i)
{
  return elements[i];
}

bool Migration::has(MeshEntity* e)
{
  return mesh->hasTag(e,tag);
}

void Migration::send(MeshEntity* e, int to)
{
  if (!has(e))
    elements.push_back(e);
  mesh->setIntTag(e,tag,&to);
}

int Migration::sending(MeshEntity* e)
{
  int to;
  mesh->getIntTag(e,tag,&to);
  return to;
}

void removeTagFromDimension(Mesh* m, MeshTag* tag, int d)
{
  MeshIterator* it = m->begin(d);
  MeshEntity* e;
  while ((e = m->iterate(it)))
    if (m->hasTag(e,tag))
      m->removeTag(e,tag);
  m->end(it);
}

int findIn(MeshEntity** a, int n, MeshEntity* e)
{
  for (int i=0; i < n; ++i)
    if (a[i]==e)
      return i;
  return -1;
}

/* returns true if the arrays have the same entities,
   regardless of ordering */
static bool sameContent(int n, MeshEntity** a, MeshEntity** b)
{
  for (int i=0; i < n; ++i)
    if (findIn(b,n,a[i])<0)
      return false;
  return true;
}

MeshEntity* findUpward(Mesh* m, int type, MeshEntity** down)
{
  if ( ! down[0]) return 0;
  Up ups;
  m->getUp(down[0],ups);
  int d = Mesh::typeDimension[type];
  int nd = Mesh::adjacentCount[type][d-1];
  Downward down2;
  for (int i=0; i < ups.n; ++i)
  {
    MeshEntity* up = ups.e[i];
    if (m->getType(up)!=type) continue;
    m->getDownward(up,d-1,down2);
    if (sameContent(nd,down,down2))
      return up;
  }
  return 0;
}

static void runEdgeDown(
    ElementVertOp*,
    MeshEntity** verts,
    MeshEntity** down)
{
  down[0] = verts[0];
  down[1] = verts[1];
}

static void runTriDown(
    ElementVertOp* o,
    MeshEntity** verts,
    MeshEntity** down)
{
  Downward ev;
  for (int i=0; i < 3; ++i)
  {
    ev[0] = verts[tri_edge_verts[i][0]]; ev[1] = verts[tri_edge_verts[i][1]];
    down[i] = o->run(Mesh::EDGE,ev);
  }
}

static void runQuadDown(
    ElementVertOp* o,
    MeshEntity** verts,
    MeshEntity** down)
{
  Downward ev;
  for (int i=0; i < 4; ++i)
  {
    ev[0] = verts[quad_edge_verts[i][0]]; ev[1] = verts[quad_edge_verts[i][1]];
    down[i] = o->run(Mesh::EDGE,ev);
  }
}

static void runTetDown(
    ElementVertOp* o,
    MeshEntity** verts,
    MeshEntity** down)
{
  Downward fv;
  for (int i=0; i < 4; ++i)
  {
    for (int j=0; j < 3; ++j)
      fv[j] = verts[tet_tri_verts[i][j]];
    down[i] = o->run(Mesh::TRIANGLE,fv);
  }
}

static void runPrismDown(
    ElementVertOp* o,
    MeshEntity** verts,
    MeshEntity** down)
{
  Downward fv;
  for (int j=0; j < 3; ++j)
    fv[j] = verts[prism_tri_verts[0][j]];
  down[0] = o->run(Mesh::TRIANGLE,fv);
  for (int i=0; i < 3; ++i)
  {
    for (int j=0; j < 4; ++j)
      fv[j] = verts[prism_quad_verts[i][j]];
    down[i+1] = o->run(Mesh::QUAD,fv);
  }
  for (int j=0; j < 3; ++j)
    fv[j] = verts[prism_tri_verts[1][j]];
  down[4] = o->run(Mesh::TRIANGLE,fv);
}

static void runPyramidDown(
    ElementVertOp* o,
    MeshEntity** verts,
    MeshEntity** down)
{
/* first 4 verts are in order for the quad */
  down[0] = o->run(Mesh::QUAD,verts);
  for (int i=0; i < 4; ++i)
  {
    MeshEntity* tv[3];
    for (int j=0; j < 3; ++j)
      tv[j] = verts[pyramid_tri_verts[i][j]];
    down[i+1] = o->run(Mesh::TRIANGLE,tv);
  }
}

void ElementVertOp::runDown(
    int type,
    MeshEntity** verts,
    MeshEntity** down)
{
  typedef void (*RunDownFunction)(
      ElementVertOp* o,
      MeshEntity** verts,
      MeshEntity** down);
  static RunDownFunction table[Mesh::TYPES] =
  {0,//vertex
   runEdgeDown,
   runTriDown,
   runQuadDown,
   runTetDown,
   0,//hex
   runPrismDown,
   runPyramidDown
  };
  RunDownFunction runDownFunction = table[type];
  runDownFunction(this,verts,down);
}

MeshEntity* ElementVertOp::run(int type, MeshEntity** verts)
{
  Downward down;
  this->runDown(type,verts,down);
  return this->apply(type,down);
}

class ElementFinder : public ElementVertOp
{
  public:
    ElementFinder(Mesh* m)
    {
      mesh = m;
    }
    virtual MeshEntity* apply(int type, MeshEntity** down)
    {
      return findUpward(mesh,type,down);
    }
  private:
    Mesh* mesh;
};

void findTriDown(
    Mesh* m,
    MeshEntity** verts,
    MeshEntity** down)
{
  ElementFinder f(m);
  return f.runDown(Mesh::TRIANGLE,verts,down);
}

MeshEntity* findElement(
    Mesh* m,
    int type,
    MeshEntity** verts)
{
  ElementFinder f(m);
  return f.run(type,verts);
}

MeshEntity* getEdgeVertOppositeVert(Mesh* m, MeshEntity* edge, MeshEntity* v)
{
  MeshEntity* ev[2];
  m->getDownward(edge,0,ev);
  if (ev[0]==v)
    return ev[1];
  return ev[0];
}

int countEntitiesOfType(Mesh* m, int type)
{
  MeshIterator* it = m->begin(Mesh::typeDimension[type]);
  MeshEntity* e;
  int count = 0;
  while ((e = m->iterate(it)))
    if (m->getType(e)==type)
      ++count;
  m->end(it);
  return count;
}

void changeMeshShape(Mesh* m, FieldShape* newShape, bool project)
{
  Field* oldCoordinateField = m->getCoordinateField();
  VectorField* newCoordinateField = new VectorField();
  newCoordinateField->init(oldCoordinateField->getName(), m, newShape,
      new TagDataOf<double>());
  if (project)
    newCoordinateField->project(oldCoordinateField);
  m->changeCoordinateField(newCoordinateField);
}

void unfreezeFields(Mesh* m) {
  Field* f;
  for (int i=0; i<m->countFields(); i++) {
    f = m->getField(i);
    if (isFrozen(f)) {
      unfreeze(f);
    }
  }
  m->hasFrozenFields = false;
}

Copy getOtherCopy(Mesh* m, MeshEntity* s)
{
  Copies remotes;
  m->getRemotes(s,remotes);
  assert(remotes.size()==1);
  Copies::iterator it = remotes.begin();
  return Copy(it->first, it->second);
}

int getDimension(Mesh* m, MeshEntity* e)
{
  return Mesh::typeDimension[m->getType(e)];
}

bool isSimplex(int type)
{
  bool const table[Mesh::TYPES] =
  {true,//VERT
   true,//EDGE
   true,//TRI
   false,//QUAD
   true,//TET
   false,//HEX
   false,//PRISM
   false,//PYRAMID
  };
  return table[type];
}

Vector3 getLinearCentroid(Mesh* m, MeshEntity* e)
{
  apf::Downward v;
  int nv = m->getDownward(e, 0, v);
  apf::Vector3 c(0,0,0);
  apf::Vector3 tmp;
  for (int i = 0; i < nv; ++i) {
    m->getPoint(v[i], 0, tmp);
    c = c + tmp;
  }
  return c / nv;
}

int countEntitiesOn(Mesh* m, ModelEntity* me, int dim)
{
  MeshIterator* it = m->begin(dim);
  MeshEntity* e;
  int n = 0;
  while ((e = m->iterate(it)))
    if (m->toModel(e) == me)
      ++n;
  m->end(it);
  return n;
}

int countOwned(Mesh* m, int dim)
{
  MeshIterator* it = m->begin(dim);
  MeshEntity* e;
  int n = 0;
  while ((e = m->iterate(it)))
    if (m->isOwned(e))
      ++n;
  m->end(it);
  return n;
}

void printStats(Mesh* m)
{
  long n[4];
  for (int i = 0; i < 4; ++i)
    n[i] = countOwned(m, i);
  PCU_Add_Longs(n, 4);
  if (!PCU_Comm_Self())
    printf("mesh entity counts: v %ld e %ld f %ld r %ld\n",
        n[0], n[1], n[2], n[3]);
}

void warnAboutEmptyParts(Mesh* m)
{
  int emptyParts = 0;
  if (!m->count(m->getDimension()))
    ++emptyParts;
  PCU_Add_Ints(&emptyParts, 1);
  if (emptyParts && (!PCU_Comm_Self()))
    fprintf(stderr,"APF warning: %d empty parts\n",emptyParts);
}

static void getRemotesArray(Mesh* m, MeshEntity* e, CopyArray& a)
{
  Copies remotes;
  m->getRemotes(e, remotes);
  a.setSize(remotes.size());
  size_t i = 0;
  APF_ITERATE(Copies, remotes, it) {
    a[i].peer = it->first;
    a[i].entity = it->second;
    ++i;
  }
}

struct NormalSharing : public Sharing
{
  NormalSharing(Mesh* m):mesh(m) {}
  bool isOwned(MeshEntity* e)
  {
    return mesh->isOwned(e);
  }
  virtual void getCopies(MeshEntity* e,
      CopyArray& copies)
  {
    if (!mesh->isShared(e))
      return;
    getRemotesArray(mesh, e, copies);
  }
  Mesh* mesh;
};

/* okay... so previously this used a min-rank rule
   for matched copies, but thats inconsistent with
   the min-count rule in MDS for regular copies,
   and the usual partition model ignores matched
   copies.
   for now, we'll build a full neighbor count
   system into this object just to implement
   the min-count rule for matched neighbors */
struct MatchedSharing : public Sharing
{
  MatchedSharing(Mesh* m):helper(m),mesh(m)
  {
    formCountMap();
  }
  size_t getNeighborCount(int peer)
  {
    assert(countMap.count(peer));
    return countMap[peer];
  }
  bool isLess(Copy const& a, Copy const& b)
  {
    size_t ca = this->getNeighborCount(a.peer);
    size_t cb = this->getNeighborCount(b.peer);
    if (ca != cb)
      return ca < cb;
    if (a.peer != b.peer)
      return a.peer < b.peer;
    return a.entity < b.entity;
  }
  Copy getOwner(MeshEntity* e)
  {
    Copy owner(PCU_Comm_Self(), e);
    CopyArray copies;
    this->getCopies(e, copies);
    APF_ITERATE(CopyArray, copies, cit)
      if (this->isLess(*cit, owner))
        owner = *cit;
    return owner;
  }
  bool isOwned(MeshEntity* e)
  {
    Copy owner = this->getOwner(e);
    return owner.peer == PCU_Comm_Self() && owner.entity == e;
  }
  void getCopies(MeshEntity* e,
      CopyArray& copies)
  {
    mesh->getMatches(e, copies);
    if (!copies.getSize())
      helper.getCopies(e, copies);
  }
  void getNeighbors(Parts& neighbors)
  {
    MeshIterator* it = mesh->begin(0);
    MeshEntity* v;
    while ((v = mesh->iterate(it))) {
      CopyArray copies;
      this->getCopies(v, copies);
      APF_ITERATE(CopyArray, copies, cit)
        neighbors.insert(cit->peer);
    }
    mesh->end(it);
    neighbors.erase(PCU_Comm_Self());
  }
  void formCountMap()
  {
    size_t count = mesh->count(mesh->getDimension());
    countMap[PCU_Comm_Self()] = count;
    PCU_Comm_Begin();
    Parts neighbors;
    getNeighbors(neighbors);
    APF_ITERATE(Parts, neighbors, nit)
      PCU_COMM_PACK(*nit, count);
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      size_t oc;
      PCU_COMM_UNPACK(oc);
      countMap[PCU_Comm_Sender()] = oc;
    }
  }
  NormalSharing helper;
  Mesh* mesh;
  std::map<int, size_t> countMap;
};

Sharing* getSharing(Mesh* m)
{
  if (m->hasMatching())
    return new MatchedSharing(m);
  return new NormalSharing(m);
}

static void getUpBridgeAdjacent(Mesh* m, MeshEntity* origin,
    int bridgeDimension, int targetDimension,
    std::set<MeshEntity*>& result)
{
  assert(targetDimension < bridgeDimension);
  Adjacent bridges;
  m->getAdjacent(origin, bridgeDimension, bridges);
  for (size_t i = 0; i < bridges.getSize(); ++i) {
    Downward targets;
    int nt = m->getDownward(bridges[i], targetDimension, targets);
    for (int j = 0; j < nt; ++j)
      result.insert(targets[j]);
  }
}

static void getDownBridgeAdjacent(Mesh* m, MeshEntity* origin,
    int bridgeDimension, int targetDimension,
    std::set<MeshEntity*>& result)
{
  assert(targetDimension > bridgeDimension);
  std::set<MeshEntity*> s;
  Downward bridges;
  int nbridges = m->getDownward(origin, bridgeDimension, bridges);
  for (int i = 0; i < nbridges; ++i) {
    Adjacent targets;
    m->getAdjacent(bridges[i], targetDimension, targets);
    for (size_t j = 0; j < targets.getSize(); ++j)
      result.insert(targets[j]);
  }
}

void getBridgeAdjacent(Mesh* m, MeshEntity* origin,
    int bridgeDimension, int targetDimension, Adjacent& result)
{
  assert(targetDimension != bridgeDimension);
  std::set<MeshEntity*> s;
  if (targetDimension < bridgeDimension)
    getUpBridgeAdjacent(m, origin, bridgeDimension, targetDimension, s);
  else
    getDownBridgeAdjacent(m, origin, bridgeDimension, targetDimension, s);
  s.erase(origin);
  result.setSize(s.size());
  if (s.size())
    std::copy(s.begin(), s.end(), result.begin());
}

int getFirstType(Mesh* m, int dim)
{
  MeshIterator* it = m->begin(dim);
  MeshEntity* e = m->iterate(it);
  m->end(it);
  return m->getType(e);
}

static int const* getVertIndices(int type, int subtype, int which)
{
  switch (subtype) {
    case Mesh::EDGE:
    switch (type) {
      case Mesh::TRIANGLE:
        return tri_edge_verts[which];
      case Mesh::QUAD:
        return quad_edge_verts[which];
      case Mesh::TET:
        return tet_edge_verts[which];
    };
    case Mesh::TRIANGLE:
    switch (type) {
      case Mesh::TET:
        return tet_tri_verts[which];
    };
  };
  fail("getVertIndices: types not supported\n");
}

void getAlignment(Mesh* m, MeshEntity* elem, MeshEntity* boundary,
    int& which, bool& flip, int& rotate)
{
  Downward ev;
  m->getDownward(elem, 0, ev);
  Downward eb;
  int neb = m->getDownward(elem, getDimension(m, boundary), eb);
  which = findIn(eb, neb, boundary);
  Downward bv;
  int nbv = m->getDownward(boundary, 0, bv);
  int const* vi = getVertIndices(m->getType(elem), m->getType(boundary), which);
  Downward ebv;
  for (int i = 0; i < nbv; ++i)
    ebv[i] = ev[vi[i]];
  int a = findIn(bv, nbv, ebv[0]);
  int b = findIn(bv, nbv, ebv[1]);
  if (b == (a + 1) % nbv)
    flip = false;
  else {
    flip = true;
    Downward tmp;
    for (int i = 0; i < nbv; ++i)
      tmp[nbv - i - 1] = bv[i];
    for (int i = 0; i < nbv; ++i)
      bv[i] = tmp[i];
  }
  rotate = findIn(bv, nbv, ebv[0]);
}

void packString(std::string s, int to)
{
  size_t len = s.length();
  PCU_COMM_PACK(to, len);
  PCU_Comm_Pack(to, s.c_str(), len);
}

std::string unpackString()
{
  std::string s;
  size_t len;
  PCU_COMM_UNPACK(len);
  s.resize(len);
  PCU_Comm_Unpack((void*)s.c_str(), len);
  return s;
}

void packTagInfo(Mesh* m, MeshTag* t, int to)
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

void unpackTagInfo(std::string& name, int& type, int& size)
{
  name = unpackString();
  PCU_COMM_UNPACK(type);
  PCU_COMM_UNPACK(size);
}

} //namespace apf
