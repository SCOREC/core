/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include <PCU.h>
#include "apfMDS.h"
#include "mds_apf.h"
#include "apfPM.h"
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <apfPartition.h>
#include <cstring>

extern "C" {

int const mds_apf_double = apf::Mesh::DOUBLE;
int const mds_apf_int = apf::Mesh::INT;
int const mds_apf_long = apf::Mesh::LONG;

}

namespace apf {

static MeshEntity* toEnt(mds_id id)
{
  //for pointers 0 is null, but mds_id 0 is ok and -1 is null
  return reinterpret_cast<MeshEntity*>(((char*)1) + id);
}

static mds_id fromEnt(MeshEntity* e)
{
  return (reinterpret_cast<char*>(e) - ((char*)1));
}

static MeshIterator* makeIter()
{
  mds_id* p = new mds_id;
  return reinterpret_cast<MeshIterator*>(p);
}

static void freeIter(MeshIterator* it)
{
  mds_id* p = reinterpret_cast<mds_id*>(it);
  delete p;
}

static void toIter(mds_id id, MeshIterator* it)
{
  *(reinterpret_cast<mds_id*>(it)) = id;
}

static mds_id fromIter(MeshIterator* it)
{
  return *(reinterpret_cast<mds_id*>(it));
}

static int mds2apf(int t_mds)
{
  static int const table[MDS_TYPES] =
  {Mesh::VERTEX
  ,Mesh::EDGE
  ,Mesh::TRIANGLE
  ,Mesh::QUAD
  ,Mesh::PRISM
  ,Mesh::PYRAMID
  ,Mesh::TET
  ,Mesh::HEX};
  return table[t_mds];

}

static int apf2mds(int t_apf)
{
  static int const table[Mesh::TYPES] =
  {MDS_VERTEX
  ,MDS_EDGE
  ,MDS_TRIANGLE
  ,MDS_QUADRILATERAL
  ,MDS_TETRAHEDRON
  ,MDS_HEXAHEDRON
  ,MDS_WEDGE
  ,MDS_PYRAMID
  };
  return table[t_apf];
}

class MeshMDS : public Mesh2
{
  public:
    MeshMDS(gmi_model* m, int d, bool isMatched_)
    {
      init(apf::getLagrange(1));
      mds_id cap[MDS_TYPES] = {};
      mesh = mds_apf_create(m, d, cap);
      isMatched = isMatched_;
    }
    MeshMDS(gmi_model* m, Mesh* from)
    {
      init(apf::getLagrange(1));
      mds_id cap[MDS_TYPES];
      cap[MDS_VERTEX] = from->count(0);
      cap[MDS_EDGE] = from->count(1);
      cap[MDS_TRIANGLE] = countEntitiesOfType(from,TRIANGLE);
      cap[MDS_QUADRILATERAL] = countEntitiesOfType(from,QUAD);
      cap[MDS_WEDGE] = countEntitiesOfType(from,PRISM);
      cap[MDS_PYRAMID] = countEntitiesOfType(from,PYRAMID);
      cap[MDS_TETRAHEDRON] = countEntitiesOfType(from,TET);
      cap[MDS_HEXAHEDRON] = countEntitiesOfType(from,HEX);
      int d = from->getDimension();
      mesh = mds_apf_create(m,d,cap);
      isMatched = from->hasMatching();
      apf::convert(from,this);
    }
    MeshMDS(gmi_model* m, const char* pathname)
    {
      init(apf::getLagrange(1));
      mesh = mds_read_smb(m, pathname);
      isMatched = PCU_Or(!mds_net_empty(&mesh->matches));
    }
    ~MeshMDS()
    {
    }
    int getDimension()
    {
      return mesh->mds.d;
    }
    std::size_t count(int dimension)
    {
      mds_id c = 0;
      for (int t = 0; t < MDS_TYPES; ++t)
        if (mds_dim[t] == dimension)
          c += mesh->mds.n[t];
      return c;
    }
    MeshIterator* begin(int dimension)
    {
      mds_id id = mds_begin(&(mesh->mds),dimension);
      MeshIterator* it = makeIter();
      toIter(id,it);
      return it;
    }
    MeshEntity* iterate(MeshIterator* it)
    {
      mds_id id = fromIter(it);
      if (id == MDS_NONE)
        return 0;
      MeshEntity* e = toEnt(id);
      id = mds_next(&(mesh->mds),id);
      toIter(id,it);
      return e;
    }
    void end(MeshIterator* it)
    {
      freeIter(it);
    }
    bool isShared(MeshEntity* e)
    {
      return mds_get_copies(&mesh->remotes, fromEnt(e));
    }
    bool isOwned(MeshEntity* e)
    {
      return getId() == getOwner(e);
    }
    int getOwner(MeshEntity* e)
    {
      void* vp = mds_get_part(mesh, fromEnt(e));
      PME* p = static_cast<PME*>(vp);
      return p->owner;
    }
    void getAdjacent(MeshEntity* e, int dimension, Adjacent& adjacent)
    {
      mds_set s;
      mds_id id = fromEnt(e);
      mds_get_adjacent(&(mesh->mds),id,dimension,&s);
      adjacent.setSize(s.n);
      for (int i = 0; i < s.n; ++i)
        adjacent[i] = toEnt(s.e[i]);
    }
    int getDownward(MeshEntity* e, int dimension, MeshEntity** adjacent)
    {
      assert((0 <= dimension) && (dimension <= 3));
      mds_set s;
      mds_id id = fromEnt(e);
      mds_get_adjacent(&(mesh->mds),id,dimension,&s);
      for (int i = 0; i < s.n; ++i)
        adjacent[i] = toEnt(s.e[i]);
      return s.n;
    }
    int countUpward(MeshEntity* e)
    {
      mds_set s;
      mds_id id = fromEnt(e);
      mds_get_adjacent(&(mesh->mds),id,mds_dim[mds_type(id)] + 1,&s);
      return s.n;
    }
    MeshEntity* getUpward(MeshEntity* e, int i)
    {
      mds_set s;
      mds_id id = fromEnt(e);
      mds_get_adjacent(&(mesh->mds),id,mds_dim[mds_type(id)] + 1,&s);
      assert(i < s.n);
      return toEnt(s.e[i]);
    }
    void getUp(MeshEntity* e, Up& up)
    {
      mds_set s;
      mds_id id = fromEnt(e);
      mds_get_adjacent(&(mesh->mds),id,mds_dim[mds_type(id)] + 1,&s);
      up.n = s.n;
      for (int i = 0; i < s.n; ++i)
        up.e[i] = toEnt(s.e[i]);
    }
    bool hasUp(MeshEntity* e)
    {
      return mds_has_up(&(mesh->mds),fromEnt(e));
    }
    void getPoint_(MeshEntity* e, int, Vector3& point)
    {
      mds_id id = fromEnt(e);
      point = Vector3(mds_apf_point(mesh,id));
    }
    void setPoint_(MeshEntity* e, int, Vector3 const& p)
    {
      mds_id id = fromEnt(e);
      p.toArray(mds_apf_point(mesh,id));
    }
    void getParam(MeshEntity* e, Vector3& p)
    {
      mds_id id = fromEnt(e);
      double* p2 = mds_apf_param(mesh,id);
      p[0] = p2[0];
      p[1] = p2[1];
    }
    void setParam(MeshEntity* e, Vector3 const& p)
    {
      mds_id id = fromEnt(e);
      double* p2 = mds_apf_param(mesh,id);
      p2[0] = p[0];
      p2[1] = p[1];
    }
    int getType(MeshEntity* e)
    {
      return mds2apf(mds_type(fromEnt(e)));
    }
    void getRemotes(MeshEntity* e, Copies& remotes)
    {
      if (!isShared(e))
        return;
      mds_copies* c = mds_get_copies(&mesh->remotes, fromEnt(e));
      for (int i = 0; i < c->n; ++i)
        remotes[c->c[i].p] = toEnt(c->c[i].e);
    }
    void getResidence(MeshEntity* e, Parts& residence)
    {
      void* vp = mds_get_part(mesh, fromEnt(e));
      PME* p = static_cast<PME*>(vp);
      for (size_t i = 0; i < p->ids.size(); ++i)
        residence.insert(p->ids[i]);
    }
    MeshTag* createDoubleTag(const char* name, int size)
    {
      mds_tag* tag;
      assert(!mds_find_tag(&mesh->tags, name));
      tag = mds_create_tag(&(mesh->tags),name,
          sizeof(double)*size, Mesh::DOUBLE);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* createIntTag(const char* name, int size)
    {
      mds_tag* tag;
      assert(!mds_find_tag(&mesh->tags, name));
      tag = mds_create_tag(&(mesh->tags),name,
          sizeof(int)*size, Mesh::INT);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* createLongTag(const char* name, int size)
    {
      mds_tag* tag;
      assert(!mds_find_tag(&mesh->tags, name));
      tag = mds_create_tag(&(mesh->tags),name,
          sizeof(long)*size, Mesh::LONG);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* findTag(const char* name)
    {
      mds_tag* tag;
      tag = mds_find_tag(&(mesh->tags),name);
      return reinterpret_cast<MeshTag*>(tag);
    }
    void destroyTag(MeshTag* t)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      mds_destroy_tag(&(mesh->tags),tag);
    }
    void getTags(DynamicArray<MeshTag*>& tags)
    {
      int c = 0;
      for (mds_tag* t = mesh->tags.first; t; t = t->next)
        ++c;
      tags.setSize(c);
      c = 0;
      for (mds_tag* t = mesh->tags.first; t; t = t->next)
        tags[c++] = reinterpret_cast<MeshTag*>(t);
    }
    void getTag(MeshEntity* e, MeshTag* t, void* data)
    {
      if (!hasTag(e,t)) {
        fprintf(stderr, "expected tag \"%s\" on entity type %d\n",
            getTagName(t), getType(e));
        abort();
      }
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      mds_id id = fromEnt(e);
      memcpy(data,mds_get_tag(tag,id),tag->bytes);
    }
    void setTag(MeshEntity* e, MeshTag* t, void const* data)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      mds_id id = fromEnt(e);
      if ( ! mds_has_tag(tag,id))
        mds_give_tag(tag,&(mesh->mds),id);
      memcpy(mds_get_tag(tag,id),data,tag->bytes);
    }
    void getDoubleTag(MeshEntity* e, MeshTag* tag, double* data)
    {
      getTag(e,tag,data);
    }
    void setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data)
    {
      setTag(e,tag,data);
    }
    void getIntTag(MeshEntity* e, MeshTag* tag, int* data)
    {
      getTag(e,tag,data);
    }
    void setIntTag(MeshEntity* e, MeshTag* tag, int const* data)
    {
      setTag(e,tag,data);
    }
    void getLongTag(MeshEntity* e, MeshTag* tag, long* data)
    {
      getTag(e,tag,data);
    }
    void setLongTag(MeshEntity* e, MeshTag* tag, long const* data)
    {
      setTag(e,tag,data);
    }
    void removeTag(MeshEntity* e, MeshTag* t)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      mds_id id = fromEnt(e);
      mds_take_tag(tag,id);
    }
    bool hasTag(MeshEntity* e, MeshTag* t)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      mds_id id = fromEnt(e);
      return mds_has_tag(tag,id);
    }
    int getTagType(MeshTag* t)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      return tag->user_type;
    }
    int getTagSize(MeshTag* t)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      size_t const table[3] =
      {sizeof(double),sizeof(int),sizeof(long)};
      return tag->bytes / table[tag->user_type];
    }
    const char* getTagName(MeshTag* t)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      return tag->name;
    }
    ModelEntity* toModel(MeshEntity* e)
    {
      return reinterpret_cast<ModelEntity*>(mds_apf_model(mesh, fromEnt(e)));
    }
    gmi_model* getModel()
    {
      return mesh->user_model;
    }
    void acceptChanges()
    {
      updateOwners(this, parts);
    }
    void migrate(Migration* plan)
    {
      apf::migrate(this,plan);
    }
    int getId()
    {
      return PCU_Comm_Self();
    }
    void writeNative(const char* fileName)
    {
      double t0 = PCU_Time();
      mesh = mds_write_smb(mesh, fileName);
      double t1 = PCU_Time();
      if (!PCU_Comm_Self())
        printf("mesh %s written in %f seconds\n", fileName, t1 - t0);
    }
    void destroyNative()
    {
      while (this->countFields())
        apf::destroyField(this->getField(0));
      while (this->countNumberings())
        apf::destroyNumbering(this->getNumbering(0));
      changeCoordinateField(0);
      gmi_model* model = static_cast<gmi_model*>(mesh->user_model);
      PCU_Thrd_Barrier();
      if (!PCU_Thrd_Self())
        gmi_destroy(model);
      PCU_Thrd_Barrier();
      mds_apf_destroy(mesh);
      mesh = 0;
    }
    void verify()
    {
      apf::verify(this);
    }
    void setRemotes(MeshEntity* e, Copies& remotes)
    {
      mds_id id = fromEnt(e);
      if (!remotes.size())
        return mds_set_copies(&mesh->remotes, &mesh->mds, id, NULL);
      mds_copies* c = mds_make_copies(remotes.size());
      c->n = 0;
      APF_ITERATE(Copies, remotes, it) {
        c->c[c->n].p = it->first;
        c->c[c->n].e = fromEnt(it->second);
        ++c->n;
      }
      mds_set_copies(&mesh->remotes, &mesh->mds, id, c);
    }
    void addRemote(MeshEntity* e, int p, MeshEntity* r)
    {
      mds_copy c;
      c.e = fromEnt(r);
      c.p = p;
      mds_add_copy(&mesh->remotes, &mesh->mds, fromEnt(e), c);
    }
    void setResidence(MeshEntity* e, Parts& residence)
    {
      mds_id id = fromEnt(e);
      PME* p = getPME(parts, residence);
      void* vp = static_cast<void*>(p);
      void* ovp = mds_get_part(mesh, id);
      if (ovp) { /* partition model classification can be NULL during
        early mesh initialization, such as after reading SMB or in
        createEntity_ */
        PME* op = static_cast<PME*>(ovp);
        putPME(parts, op);
      }
      mds_set_part(mesh, id, vp);
    }
    void increment(MeshIterator* it)
    {
      toIter(mds_next(&(mesh->mds),fromIter(it)),it);
    }
    bool isDone(MeshIterator* it)
    {
      return fromIter(it) == MDS_NONE;
    }
    MeshEntity* deref(MeshIterator* it)
    {
      return toEnt(fromIter(it));
    }
    MeshEntity* createVert_(ModelEntity* c)
    {
      return createEntity_(VERTEX,c,0);
    }
    MeshEntity* createEntity_(int type, ModelEntity* c,
                                      MeshEntity** down)
    {
      int t = apf2mds(type);
      int dim = mds_dim[t];
      if (dim > mesh->mds.d) {
        fprintf(stderr,"error: creating entity of dimension %d "
                       "in mesh of dimension %d\n", dim, mesh->mds.d);
        fprintf(stderr,"please use apf::changeMdsDimension\n");
        abort();
      }
      mds_set s;
      if (type != VERTEX) {
        s.n = mds_degree[t][mds_dim[t]-1];
        for (int i = 0; i < s.n; ++i)
          s.e[i] = fromEnt(down[i]);
      }
      mds_id id = mds_apf_create_entity(
          mesh, t, reinterpret_cast<gmi_ent*>(c), s.e);
      MeshEntity* e = toEnt(id);
      apf::Parts r;
      r.insert(getId());
      setResidence(e, r);
      return e;
    }
    void destroy_(MeshEntity* e)
    {
      mds_id id = fromEnt(e);
      void* ovp = mds_get_part(mesh, id);
      PME* op = static_cast<PME*>(ovp);
      putPME(parts, op);
      mds_apf_destroy_entity(mesh,id);
    }
    bool hasMatching()
    {
      return isMatched;
    }
    void getMatches(MeshEntity* e, Matches& m)
    {
      mds_copies* c = mds_get_copies(&mesh->matches, fromEnt(e));
      if (!c)
        return;
      m.setSize(c->n);
      for (int i = 0; i < c->n; ++i) {
        m[i].entity = toEnt(c->c[i].e);
        m[i].peer = c->c[i].p;
      }
    }
    void addMatch(MeshEntity* e, int peer, MeshEntity* match)
    {
      mds_copy c;
      c.e = fromEnt(match);
      c.p = peer;
      mds_add_copy(&mesh->matches, &mesh->mds, fromEnt(e), c);
    }
    void clearMatches(MeshEntity* e)
    {
      mds_set_copies(&mesh->matches, &mesh->mds, fromEnt(e), 0);
    }
    double getElementBytes(int type)
    {
      static double const table[TYPES] =
      {1  , //vertex
       1  , //edge
       1  , //triangle
       1  , //quad
       250, //tet
       1  , //hex
       350, //prism
       300, //pyramid
      };
      return table[type];
    }
    mds_apf* mesh;
    PM parts;
    bool isMatched;
};

Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched)
{
  Mesh2* m = new MeshMDS(model, dim, isMatched);
  initResidence(m, dim);
  return m;
}

Mesh2* createMdsMesh(gmi_model* model, Mesh* from)
{
  return new MeshMDS(model, from);
}

Mesh2* loadMdsMesh(gmi_model* model, const char* meshfile)
{
  double t0 = PCU_Time();
  Mesh2* m = new MeshMDS(model, meshfile);
  initResidence(m, m->getDimension());
  stitchMesh(m);
  m->acceptChanges();
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("mesh %s loaded in %f seconds\n", meshfile, t1 - t0);
  printStats(m);
  warnAboutEmptyParts(m);
  return m;
}

Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile)
{
  double t0 = PCU_Time();
  PCU_Thrd_Barrier();
  static gmi_model* model;
  if (!PCU_Thrd_Self())
    model = gmi_load(modelfile);
  PCU_Thrd_Barrier();
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("model %s loaded in %f seconds\n", modelfile, t1 - t0);
  return loadMdsMesh(model, meshfile);
}

void reorderMdsMesh(Mesh2* mesh)
{
  MeshMDS* m = static_cast<MeshMDS*>(mesh);
  m->mesh = mds_reorder(m->mesh);
}

static int globalFactor;
static Mesh2* globalMesh;
static Migration* globalPlan;
static void (*globalThrdCall)(Mesh2*);

Mesh2* repeatMdsMesh(Mesh2* m, gmi_model* g, Migration* plan, int factor)
{
  double t0 = PCU_Time();
  int self = PCU_Comm_Self();
  bool isOriginal = ((PCU_Comm_Self() % factor) == 0);
  int dim;
  bool isMatched;
  PCU_Comm_Begin();
  if (isOriginal) {
    dim = m->getDimension();
    isMatched = m->hasMatching();
    for (int i = 1; i < factor; ++i) {
      int clone = self + i;
      PCU_COMM_PACK(clone, dim);
      PCU_COMM_PACK(clone, isMatched);
      packDataClone(m, clone);
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(dim);
    PCU_COMM_UNPACK(isMatched);
    m = makeEmptyMdsMesh(g, dim, isMatched);
    unpackDataClone(m);
  }
  apf::Multiply remap(factor);
  apf::remapPartition(m, remap);
  if (!isOriginal)
    plan = new apf::Migration(m, m->findTag("apf_migrate"));
  m->migrate(plan);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("mesh repartitioned from %d to %d in %f seconds\n",
        PCU_Comm_Peers() / factor,
        PCU_Comm_Peers(),
        t1 - t0);
  return m;
}

extern "C" void* splitThrdMain(void*)
{
  Mesh2* m;
  m = repeatMdsMesh(globalMesh, globalMesh->getModel(), globalPlan, globalFactor);
  globalThrdCall(m);
  return NULL;
}

void splitMdsMesh(Mesh2* m, Migration* plan, int n, void (*runAfter)(Mesh2*))
{
  globalFactor = n;
  globalMesh = m;
  globalPlan = plan;
  globalThrdCall = runAfter;
  PCU_Thrd_Run(n, splitThrdMain, NULL);
}

bool alignMdsMatches(Mesh2* in)
{
  if (!in->hasMatching())
    return false;
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_align_matches(m->mesh);
}

bool alignMdsRemotes(Mesh2* in)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_align_remotes(m->mesh);
}

void deriveMdsModel(Mesh2* in)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_derive_model(m->mesh);
}

void changeMdsDimension(Mesh2* in, int d)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_change_dimension(&(m->mesh->mds), d);
}

int getMdsIndex(Mesh2* in, MeshEntity* e)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  mds* mds = &(m->mesh->mds);
  int i = 0;
  mds_id id = fromEnt(e);
  int type = mds_type(id);
  for (int t = 0; t < type; ++t)
    if (mds_dim[t] == mds_dim[type])
      i += mds->n[t];
  i += mds_index(id);
  return i;
}

MeshEntity* getMdsEntity(Mesh2* in, int dimension, int index)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  mds* mds = &(m->mesh->mds);
  for (int t = 0; t < MDS_TYPES; ++t)
    if (mds_dim[t] == dimension) {
      if (index < mds->n[t])
        return toEnt(mds_identify(t, index));
      else
        index -= mds->n[t];
    }
  return 0;
}

}
