/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "apfMDS.h"
#include "mds_apf.h"
#include "apfPM.h"
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <PCU.h>
#include <cstring>

extern "C" {

int const mds_apf_double = apf::Mesh::DOUBLE;
int const mds_apf_int = apf::Mesh::INT;

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

int getMdsId(MeshEntity* e)
{
  return mds_index(fromEnt(e));
}

static MeshIterator* makeIter()
{
  return static_cast<MeshIterator*>(malloc(sizeof(mds_id)));
}

static void freeIter(MeshIterator* it)
{
  free(static_cast<void*>(it));
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
  ,Mesh::TET};
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
  ,-1
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
      if (!isShared(e))
        return true;
      mds_copies* c = mds_get_copies(&mesh->remotes, fromEnt(e));
      return getId() < c->c[0].p;
    }
    int getOwner(MeshEntity* e)
    {
      if (isOwned(e))
        return getId();
      return mds_get_copies(&mesh->remotes, fromEnt(e))->c[0].p;
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
    void getPoint_(MeshEntity* e, int node, Vector3& point)
    {
      mds_id id = fromEnt(e);
      point = Vector3(mds_apf_point(mesh,id));
    }
    void setPoint_(MeshEntity* e, int node, Vector3 const& p)
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
      tag = mds_create_tag(&(mesh->tags),&(mesh->mds),name,
          sizeof(double)*size, Mesh::DOUBLE);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* createIntTag(const char* name, int size)
    {
      mds_tag* tag;
      assert(!mds_find_tag(&mesh->tags, name));
      tag = mds_create_tag(&(mesh->tags),&(mesh->mds),name,
          sizeof(int)*size, Mesh::INT);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* createLongTag(const char* name, int size)
    {
      mds_tag* tag;
      assert(!mds_find_tag(&mesh->tags, name));
      tag = mds_create_tag(&(mesh->tags),&(mesh->mds),name,
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
      assert(hasTag(e,t));
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
    int getModelType(ModelEntity* e)
    {
      return mds_model_dim(mesh, reinterpret_cast<gmi_ent*>(e));
    }
    int getModelTag(ModelEntity* e)
    {
      return mds_model_id(mesh, reinterpret_cast<gmi_ent*>(e));
    }
    ModelEntity* findModelEntity(int type, int tag)
    {
      return reinterpret_cast<ModelEntity*>(mds_find_model(mesh,type,tag));
    }
    ModelEntity* toModel(MeshEntity* e)
    {
      return reinterpret_cast<ModelEntity*>(mds_apf_model(mesh,fromEnt(e)));
    }
    void getModelFaceNormal(ModelEntity* face, Vector3 const& p,
                                    Vector3& n)
    {
      abort();
    }
    void getModelEdgeTangent(ModelEntity* edge, double p, Vector3& n)
    {
      abort();
    }
    bool snapToModel(ModelEntity* m, Vector3 const& p, Vector3& x)
    {
      gmi_eval(mesh->user_model,
               reinterpret_cast<gmi_ent*>(m),
               &p[0], &x[0]);
      return true;
    }
    void preMigrate_()
    {
    }
    void postMigrate_()
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
      mesh = mds_write_smb(mesh, fileName);
    }
    void destroyNative()
    {
      while (this->countFields())
        apf::destroyField(this->getField(0));
      while (this->countNumberings())
        apf::destroyNumbering(this->getNumbering(0));
      changeCoordinateField(0);
      gmi_model* model = static_cast<gmi_model*>(mesh->user_model);
      PCU_Barrier();
      if (!PCU_Thrd_Self())
        gmi_destroy(model);
      PCU_Barrier();
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
      PME* p = getPME(parts, residence);
      void* vp = static_cast<void*>(p);
      mds_set_part(mesh, fromEnt(e), vp);
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
    bool canSnap()
    {
      return gmi_can_eval(mesh->user_model);
    }
    void getParamOn(ModelEntity* g, MeshEntity* e, Vector3& p)
    {
      ModelEntity* from_g = toModel(e);
      if (g == from_g)
        return getParam(e, p);
      gmi_ent* from = reinterpret_cast<gmi_ent*>(from_g);
      gmi_ent* to = reinterpret_cast<gmi_ent*>(g);
      Vector3 from_p;
      getParam(e, from_p);
      gmi_reparam(mesh->user_model, from, &from_p[0], to, &p[0]);
    }
    bool getPeriodicRange(ModelEntity* g, int axis, double range[2])
    {
      gmi_ent* e = reinterpret_cast<gmi_ent*>(g);
      gmi_range(mesh->user_model, e, axis, range);
      return gmi_periodic(mesh->user_model, e, axis);
    }
    MeshEntity* createVert_(ModelEntity* c)
    {
      return createEntity_(VERTEX,c,0);
    }
    MeshEntity* createEntity_(int type, ModelEntity* c,
                                      MeshEntity** down)
    {
      int t = apf2mds(type);
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
      setResidence(e,r);
      return e;
    }
    void destroy_(MeshEntity* e)
    {
      mds_id id = fromEnt(e);
      mds_apf_destroy_entity(mesh,id);
    }
    void stitch()
    {
      stitchMesh(this);
      updateOwners(this, parts);
    }
    bool hasMatching()
    {
      return isMatched;
    }
    void getMatches(MeshEntity* e, Matches& m)
    {
    }
    void addMatch(MeshEntity* e, int peer, MeshEntity* match)
    {
      abort();
    }
    void clearMatches(MeshEntity* e)
    {
    }
    void repartition(MeshTag* elementWeights, double maximumImbalance)
    {
      abort();
    }
    gmi_model* getPumiModel()
    {
      return static_cast<gmi_model*>(mesh->user_model);
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

Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile)
{
  PCU_Barrier();
  static gmi_model* model;
  if (!PCU_Thrd_Self())
    model = gmi_load(modelfile);
  PCU_Barrier();
  Mesh2* m = new MeshMDS(model, meshfile);
  initResidence(m, m->getDimension());
  m->stitch();
  return m;
}

void defragMdsMesh(Mesh2* mesh)
{
  MeshMDS* m = static_cast<MeshMDS*>(mesh);
  m->mesh = mds_reorder(m->mesh);
}

gmi_model* getMdsModel(Mesh2* mesh)
{
  return static_cast<MeshMDS*>(mesh)->getPumiModel();
}

}
