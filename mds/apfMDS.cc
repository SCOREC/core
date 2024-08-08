/******************************************************************************

  Copyright 2014 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include <PCUObj.h>
#include <lionPrint.h>
#include "apfMDS.h"
#include "mds_apf.h"
#include "apfPM.h"
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <apfPartition.h>
#include <apfFile.h>
#include <cstring>
#include <pcu_util.h>
#include <cstdlib>
#include <stdint.h>
#include <limits>
#include <deque>

extern "C" {

int const mds_apf_double = apf::Mesh::DOUBLE;
int const mds_apf_int = apf::Mesh::INT;
int const mds_apf_long = apf::Mesh::LONG;

}

typedef std::map<std::pair<int, int>, int> ModelEdgeTags;

namespace apf {

static int getModelEdgeTag(int facetag1, int facetag2, ModelEdgeTags &tags,
                    long *minAvbl)
{
  std::pair<int, int> key;
  key.first = (facetag1 < facetag2) ? facetag1 : facetag2;
  key.second = (facetag1 < facetag2) ? facetag2 : facetag1;
  if (tags.count(key) == 0) {
    tags[key] = *minAvbl;
    (*minAvbl)++;
  }
  return tags[key];
}

static int getFaceIdInRegion(apf::Mesh* mesh, apf::MeshEntity* region,
                      int* bface_data)
{
  apf::Downward verts;
  apf::MeshTag* vIDTag = mesh->findTag("_vert_id");
  int vID;
  mesh->getDownward(region, 0, verts);
  // Go through all vertices. What vertex is not on the face can be used to determine the face id.
  // TODO: Good way to assert that the rest of the 3 actually exist?
  mesh->getIntTag(verts[0], vIDTag, &vID);
  if (vID != bface_data[2] && vID != bface_data[3] && vID != bface_data[4])
    return 2;
  mesh->getIntTag(verts[1], vIDTag, &vID);
  if (vID != bface_data[2] && vID != bface_data[3] && vID != bface_data[4])
    return 3;
  mesh->getIntTag(verts[2], vIDTag, &vID);
  if (vID != bface_data[2] && vID != bface_data[3] && vID != bface_data[4])
    return 1;
  mesh->getIntTag(verts[3], vIDTag, &vID);
  if (vID != bface_data[2] && vID != bface_data[3] && vID != bface_data[4])
    return 0;
  return 12; // Should give segmentation fault
}

static int getEdgeIdInFace(apf::Mesh* mesh, apf::MeshEntity* face,
                      int* bedge_data)
{
  apf::Downward verts, edges;
  apf::MeshTag* vIDTag = mesh->findTag("_vert_id");
  int vID[2], eID;
  mesh->getDownward(face, 1, edges);
  for (eID = 0; eID < 3; ++eID) {
    mesh->getDownward(edges[eID], 0, verts);
    mesh->getIntTag(verts[0], vIDTag, &vID[0]);
    mesh->getIntTag(verts[1], vIDTag, &vID[1]);
    if((vID[0] == bedge_data[2] && vID[1] == bedge_data[3]) ||
       (vID[0] == bedge_data[3] && vID[1] == bedge_data[2])) {
      return eID;
    }
  }

  return 12; // Should give segmentation fault
}

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

static Mesh::Type mds2apf(int t_mds)
{
  static Mesh::Type const table[MDS_TYPES] =
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
    MeshMDS()
    {
      mesh = 0;
      isMatched = false;
      ownsModel = false;
    }
    MeshMDS(gmi_model* m, int d, bool isMatched_, pcu::PCU *PCUObj)
    {
      init(apf::getLagrange(1), PCUObj);
      mds_id cap[MDS_TYPES] = {};
      mesh = mds_apf_create(m, d, cap);
      isMatched = isMatched_;
      ownsModel = true;
    }
    MeshMDS(gmi_model* m, Mesh* from, 
            apf::MeshEntity** nodes, apf::MeshEntity** elems, bool copy_data=true)
    {
      init(apf::getLagrange(1), from->getPCU());
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
      ownsModel = true;
      apf::convert(from, this, nodes, elems, copy_data);
    }

    MeshMDS(gmi_model* m, const char* pathname, pcu::PCU *PCUObj)
    {
      init(apf::getLagrange(1), PCUObj);
      mesh = mds_read_smb(this->getPCU()->GetCHandle(), m, pathname, 0, this);
      isMatched = this->getPCU()->Or(!mds_net_empty(&mesh->matches));
      ownsModel = true;
    }
    ~MeshMDS()
    {
      if (mesh)
        destroyNative();
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

    // return true if adjacency *from_dim <--> to_dim*  is stored
    bool hasAdjacency(int from_dim, int to_dim)
    {
      return (mesh->mds.mrm[from_dim][to_dim] == 1);
    }

    // store adjacency *from_dim <--> to_dim* if not stored
    void createAdjacency(int from_dim, int to_dim)
    {
      if (mesh->mds.mrm[from_dim][to_dim] != 1)
        mds_add_adjacency(&(mesh->mds), from_dim, to_dim);
    }
    // remove adjacency *from_dim <--> to_dim* except for one-level apart adjacency
    void deleteAdjacency(int from_dim, int to_dim)
    {
      if (mesh->mds.mrm[from_dim][to_dim] == 1 && (abs(from_dim-to_dim)>1))
        mds_remove_adjacency(&(mesh->mds), from_dim, to_dim);
    }
    bool isShared(MeshEntity* e)
    {
      return mds_get_copies(&mesh->remotes, fromEnt(e));
    }
    bool isGhost(MeshEntity* e)
    {
      MeshTag* t = findTag("ghost_tag");
      if (t && hasTag(e, t))
        return true;
      return false;
    }

    void deleteGhost(MeshEntity* e)
    {
      mds_set_copies(&mesh->ghosts, &mesh->mds, fromEnt(e), NULL);
    }

    bool isGhosted(MeshEntity* e)
    {
      MeshTag* t = findTag("ghosted_tag");
      if (t && hasTag(e, t))
        return true;
      return false;
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
      PCU_ALWAYS_ASSERT((0 <= dimension) && (dimension <= 3));
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
      PCU_ALWAYS_ASSERT(i < s.n);
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
    Type getType(MeshEntity* e)
    {
      return mds2apf(mds_type(fromEnt(e)));
    }
    void getRemotes(MeshEntity* e, Copies& remotes)
    {
      if (!isShared(e))
        return;
      mds_copies* c = mds_get_copies(&mesh->remotes, fromEnt(e));
      PCU_ALWAYS_ASSERT(c != NULL);
      for (int i = 0; i < c->n; ++i)
        remotes[c->c[i].p] = toEnt(c->c[i].e);
    }

// seol
    int getGhosts(MeshEntity* e, Copies& ghosts)
    {
      mds_copies* c = mds_get_copies(&mesh->ghosts, fromEnt(e));
      if (c==NULL) return 0;
      for (int i = 0; i < c->n; ++i)
        ghosts[c->c[i].p] = toEnt(c->c[i].e);
      return c->n;
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
      PCU_ALWAYS_ASSERT(!mds_find_tag(&mesh->tags, name));
      tag = mds_create_tag(&(mesh->tags),name,
          sizeof(double)*size, Mesh::DOUBLE);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* createIntTag(const char* name, int size)
    {
      mds_tag* tag;
      PCU_ALWAYS_ASSERT_VERBOSE(!mds_find_tag(&mesh->tags, name), name);
      tag = mds_create_tag(&(mesh->tags),name,
          sizeof(int)*size, Mesh::INT);
      return reinterpret_cast<MeshTag*>(tag);
    }
    MeshTag* createLongTag(const char* name, int size)
    {
      mds_tag* tag;
      PCU_ALWAYS_ASSERT(!mds_find_tag(&mesh->tags, name));
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
        lion_eprint(1, "expected tag \"%s\" on entity type %d\n",
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
    void renameTag(MeshTag* t, const char* newName)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      mds_rename_tag(tag,newName);
    }
    /* \brief 16 bit additive checksum of a tag
     * \remark the code is from
     *  http://barrgroup.com/Embedded-Systems/How-To/Additive-Checksums
     */
    unsigned getTagChecksum(MeshTag* t, int type)
    {
      mds_tag* tag;
      tag = reinterpret_cast<mds_tag*>(t);
      /* count the number of 'live' indices */
      int numLive = 0;
      for (int i=0; i < mesh->mds.end[type]; ++i) {
        if (mesh->mds.free[type][i] == MDS_LIVE)
          numLive++;
      }
      int nWords = numLive / sizeof(uint16_t);
      uint16_t* data = reinterpret_cast<uint16_t*>(tag->data[type]);
      uint32_t sum = 0;
      while (nWords-- > 0)
        sum += *(data++);
      /* Use carries to compute 1's complement sum. */
      sum = (sum >> 16) + (sum & 0xFFFF);
      sum += sum >> 16;
      /* Return the inverted 16-bit result.  */
      return ((unsigned) ~sum);
    }
    ModelEntity* toModel(MeshEntity* e)
    {
      return reinterpret_cast<ModelEntity*>(mds_apf_model(mesh, fromEnt(e)));
    }
    gmi_model* getModel()
    {
      return mesh->user_model;
    }
    void setModel(gmi_model* newModel)
    {
      mesh->user_model = newModel;
      return;
    }
    void acceptChanges()
    {
      updateOwners(this, pmodel);
    }

    void migrate(Migration* plan)
    {
      apf::migrate(this,plan);
    }
    int getId()
    {
      return this->getPCU()->Self();
    }
    void writeNative(const char* fileName)
    {
      double t0 = pcu::Time();
      mesh = mds_write_smb(this->getPCU()->GetCHandle(), mesh, fileName, 0, this);
      double t1 = pcu::Time();
      if (!this->getPCU()->Self())
        lion_oprint(1,"mesh %s written in %f seconds\n", fileName, t1 - t0);
    }
    void destroyNative()
    {
      while (this->countFields())
        apf::destroyField(this->getField(0));
      while (this->countNumberings())
        apf::destroyNumbering(this->getNumbering(0));
      while (this->countGlobalNumberings())
        apf::destroyGlobalNumbering(this->getGlobalNumbering(0));
      apf::destroyField(coordinateField);
      coordinateField = 0;
      gmi_model* model = static_cast<gmi_model*>(mesh->user_model);
      if (ownsModel && model)
        gmi_destroy(model);
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

//seol
    void clearRemotes(MeshEntity* e)
    {
      mds_set_copies(&mesh->remotes, &mesh->mds, fromEnt(e), 0);
    }

    void addGhost(MeshEntity* e, int p, MeshEntity* r)
    {
      mds_copy c;
      c.e = fromEnt(r);
      c.p = p;
      mds_add_copy(&mesh->ghosts, &mesh->mds, fromEnt(e), c);
    }

    void setResidence(MeshEntity* e, Parts& residence)
    {
      mds_id id = fromEnt(e);
      PME* p = getPME(pmodel, residence);
      void* vp = static_cast<void*>(p);
      void* ovp = mds_get_part(mesh, id);
      if (ovp) { /* partition model classification can be NULL during
        early mesh initialization, such as after reading SMB or in
        createEntity_ */
        PME* op = static_cast<PME*>(ovp);
        putPME(pmodel, op);
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
        lion_eprint(1,"error: creating entity of dimension %d "
                       "in mesh of dimension %d\n", dim, mesh->mds.d);
        lion_eprint(1,"please use apf::changeMdsDimension\n");
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
      if (ovp)
      {
        PME* op = static_cast<PME*>(ovp);
        putPME(pmodel, op);
      }
      mds_apf_destroy_entity(mesh,id);
    }

    void setModelEntity(MeshEntity* e, ModelEntity* c)
    {
      mds_apf_set_model(mesh, fromEnt(e),
         reinterpret_cast<gmi_ent*>(c));
    }
    bool hasMatching()
    {
      return isMatched;
    }
    void getMatches(MeshEntity* e, Matches& m)
    {
      mds_copies* c = mds_get_copies(&mesh->matches, fromEnt(e));
      if (!c) {
        m.setSize(0);
        return;
      }
      m.setSize(c->n);
      for (int i = 0; i < c->n; ++i) {
        m[i].entity = toEnt(c->c[i].e);
        m[i].peer = c->c[i].p;
      }
    }
    void getDgCopies(MeshEntity* e, DgCopies& dgCopies, ModelEntity* me)
    {
      (void) e;
      (void) dgCopies;
      (void) me;
      PCU_ALWAYS_ASSERT_VERBOSE(false, "error: getDgCopies for MDS is not implemented yet! ");
    }
    void addMatch(MeshEntity* e, int peer, MeshEntity* match)
    {
      PCU_ALWAYS_ASSERT(isMatched);
      mds_copy c;
      c.e = fromEnt(match);
      c.p = peer;
      mds_add_copy(&mesh->matches, &mesh->mds, fromEnt(e), c);
    }
    void clearMatches(MeshEntity* e)
    {
      mds_set_copies(&mesh->matches, &mesh->mds, fromEnt(e), 0);
    }
    void clear_()
    {
      mesh = mds_apf_create(mesh->user_model, mesh->mds.d, mesh->mds.n);
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
    PM pmodel;
    bool isMatched;
    bool ownsModel;
};

Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched, pcu::PCU *PCUObj)
{
  Mesh2* m = new MeshMDS(model, dim, isMatched, PCUObj);
  initResidence(m, dim);
  return m;
}

// seol -- reorder input mesh before conversion
//         starting vtx: a vtx with min Y
struct Queue {
  bool has(apf::MeshEntity* e) { return h.count(e); }
  void push(apf::MeshEntity* e)
  {
    q.push_back(e);
    h.insert(e);
  }
  void pushVector(std::vector<apf::MeshEntity*> const& l)
  {
    for (size_t i = 0; i < l.size(); ++i)
      push(l[i]);
  }
  apf::MeshEntity* pop()
  {
    apf::MeshEntity* e;
    e = q.front();
    q.pop_front();
    h.erase(e);
    return e;
  }
  bool empty() { return q.empty(); }
  std::deque<apf::MeshEntity*> q;
  std::set<apf::MeshEntity*> h;
};

int classifDim(gmi_model* model, Mesh* m, apf::MeshEntity* e)
{
  gmi_ent* clas=(gmi_ent*)m->toModel(e);
  return gmi_dim(model, clas);
}

apf::MeshEntity* findFirst(apf::Mesh* m)
{
  apf::MeshEntity* v;
  apf::MeshEntity* best;
  apf::MeshIterator* it = m->begin(0);
  best = m->iterate(it);
  Vector3 coord;
  m->getPoint(best, 0, coord);
  double min_Y=coord[1];
  while ((v = m->iterate(it)))
  {
    m->getPoint(v, 0, coord);  
    if (min_Y > coord[1])
    {
      best = v;
      min_Y = coord[1];
    }
  }
  m->end(it);
  return best;
}

bool visited(Queue& q, apf::Numbering* nn, MeshEntity* e)
{
  return apf::isNumbered(nn, e, 0, 0) || q.has(e);
}

bool hasNode(Mesh* m, MeshEntity* e)
{
  if (m->getShape()->countNodesOn(m->getType(e))>0) 
    return true;
  return false;
}

Mesh2* createMdsMesh(gmi_model* model, Mesh* from, bool reorder, bool copy_data)
{
  if (!reorder)
    return new MeshMDS(model, from, NULL, NULL, copy_data);

  int mesh_dim = from->getDimension();

  // reorder and create mesh
  apf::Numbering* nn = apf::createNumbering(from, "node", getConstant(0), 1);
  apf::Numbering* en = apf::createNumbering(from, "elem", getConstant(mesh_dim), 1);

  Queue q;
  q.push(findFirst(from));

  // node and element number starts from 0
  int labelnode = 0;
  int labelelem = 0;

  std::vector<MeshEntity*> node_arr;
  std::vector<MeshEntity*> elem_arr;

  node_arr.resize(from->count(0)+1);
  elem_arr.resize(from->count(mesh_dim)+1);

  MeshEntity* otherVtx;
  MeshEntity* edge;
  MeshEntity* elem;

  while (!q.empty()) 
  {
    MeshEntity* vtx = q.pop();
    if (!apf::isNumbered(nn, vtx, 0, 0))
    {
      node_arr[labelnode] = vtx;
      apf::number(nn, vtx, 0, 0, labelnode);

      ++labelnode;
    }

    std::vector<MeshEntity*> entities;
    apf::Adjacent edges;
    from->getAdjacent(vtx,1, edges);
    for (size_t i = 0; i < edges.getSize(); ++i) 
    {
      edge = edges[i];
      apf::Adjacent adjacent;
      from->getAdjacent(edge, mesh_dim, adjacent);      
      for (size_t j = 0; j < adjacent.getSize(); ++j) 
      {
        elem = adjacent[j];
        if (!apf::isNumbered(en, elem, 0, 0))
        {
          elem_arr[labelelem] = elem;
          apf::number(en, elem, 0, 0, labelelem);
          ++labelelem;
        }
      }
      otherVtx = apf::getEdgeVertOppositeVert(from, edge, vtx);
      if (!visited(q, nn, otherVtx))
        entities.push_back(otherVtx);
    }
    q.pushVector(entities);
  } // while

  destroyNumbering(nn);
  destroyNumbering(en);

  return new MeshMDS(model, from, &(node_arr[0]), &(elem_arr[0]), copy_data);
}

Mesh2* loadSerialMdsMesh(gmi_model* model, const char* meshfile, pcu::PCU *PCUObj)
{
  Mesh2* m;
  m = new MeshMDS(model, meshfile, PCUObj);
  return m;
}

Mesh2* loadMdsMesh(gmi_model* model, const char* meshfile, pcu::PCU *PCUObj)
{
  double t0 = pcu::Time();
  Mesh2* m = new MeshMDS(model, meshfile, PCUObj);
  initResidence(m, m->getDimension());
  stitchMesh(m);
  m->acceptChanges();

  if (!m->getPCU()->Self())
    lion_oprint(1,"mesh %s loaded in %f seconds\n", meshfile, pcu::Time() - t0);
  printStats(m);
  warnAboutEmptyParts(m);
  return m;
}

Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile, pcu::PCU *PCUObj)
{
  double t0 = pcu::Time();
  static gmi_model* model;
  model = gmi_load(modelfile);
  if (!PCUObj->Self())
    lion_oprint(1,"model %s loaded in %f seconds\n", modelfile, pcu::Time() - t0);

  return loadMdsMesh(model, meshfile, PCUObj);
}

void reorderMdsMesh(Mesh2* mesh, MeshTag* t)
{
  double t0 = pcu::Time();
  MeshMDS* m = static_cast<MeshMDS*>(mesh);
  mds_tag* vert_nums;
  if (t) {
    PCU_ALWAYS_ASSERT(mesh->getTagType(t) == Mesh::INT);
    vert_nums = reinterpret_cast<mds_tag*>(t);
  } else {
    vert_nums = mds_number_verts_bfs(m->mesh);
  }
  m->mesh = mds_reorder(mesh->getPCU()->GetCHandle(), m->mesh, 0, vert_nums);
  if (!mesh->getPCU()->Self())
    lion_oprint(1,"mesh reordered in %f seconds\n", pcu::Time()-t0);
}


Mesh2* expandMdsMesh(Mesh2* m, gmi_model* g, int inputPartCount, pcu::PCU *expandedPCU)
{
  double t0 = pcu::Time();
  int self = expandedPCU->Self();
  int outputPartCount = expandedPCU->Peers();
  apf::Expand expand(inputPartCount, outputPartCount);
  apf::Contract contract(inputPartCount, outputPartCount);
  bool isOriginal = contract.isValid(self);
  int dim;
  bool isMatched;
  expandedPCU->Begin();
  if (isOriginal) {
    PCU_ALWAYS_ASSERT(m != 0);
    dim = m->getDimension();
    isMatched = m->hasMatching();
    for (int i = self + 1; i < outputPartCount && !contract.isValid(i); ++i) {
      expandedPCU->Pack(i, dim);
      expandedPCU->Pack(i, isMatched);
      packDataClone(m, i, expandedPCU);
    }
  }
  expandedPCU->Send();
  while (expandedPCU->Receive()) {
    expandedPCU->Unpack(dim);
    expandedPCU->Unpack(isMatched);
    m = makeEmptyMdsMesh(g, dim, isMatched, expandedPCU);
    unpackDataClone(m);
  }
  PCU_ALWAYS_ASSERT(m != 0);
  apf::remapPartition(m, expand);
  double t1 = pcu::Time();
  if (!m->getPCU()->Self())
    lion_oprint(1,"mesh expanded from %d to %d parts in %f seconds\n",
        inputPartCount, outputPartCount, t1 - t0);
  return m;
}

Mesh2* repeatMdsMesh(Mesh2* m, gmi_model* g, Migration* plan,
    int factor, pcu::PCU *PCUObj)
{
  m = expandMdsMesh(m, g, PCUObj->Peers() / factor, PCUObj);
  double t0 = pcu::Time();
  if (PCUObj->Self() % factor != 0)
    plan = new apf::Migration(m, m->findTag("apf_migrate"));
  m->migrate(plan);
  double t1 = pcu::Time();
  if (!PCUObj->Self())
    lion_oprint(1,"mesh migrated from %d to %d in %f seconds\n",
        PCUObj->Peers() / factor,
        PCUObj->Peers(),
        t1 - t0);
  return m;
}

bool alignMdsMatches(Mesh2* in)
{
  if (!in->hasMatching())
    return false;
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_align_matches(in->getPCU()->GetCHandle(), m->mesh);
}

bool alignMdsRemotes(Mesh2* in)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_align_remotes(in->getPCU()->GetCHandle(), m->mesh);
}

void deriveMdsModel(Mesh2* in)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  return mds_derive_model(m->mesh);
}

void deriveMdlFromManifold(Mesh2* mesh, bool* isModelVert,
                           int nBFaces, int (*bFaces)[5],
                           GlobalToVert &globalToVert,
                           std::map<int, apf::MeshEntity*> &globalToRegion)
{
  PCU_ALWAYS_ASSERT_VERBOSE(!mesh->findTag("_classifn_data"),
          "MeshTag name \"_classifn_data\" is used internally in this method\n");
  apf::MeshTag* classifnTag = mesh->createIntTag("_classifn_data", 2);
  int tagData[2], newTagData[2];
  long minAvbl = 1;
  // This is set by during apf::construct, and using anything else leads to an
  // additional region
  long DEFAULT_REGION_ID = 0;

  PCU_ALWAYS_ASSERT_VERBOSE(!mesh->findTag("_vert_id"),
          "MeshTag name \"_vert_id\" is used internally in this method\n");
  apf::MeshTag* vIDTag = mesh->createLongTag("_vert_id", 1);
  for (apf::GlobalToVert::iterator vit = globalToVert.begin();
       vit !=  globalToVert.end(); vit++) {
    mesh->setLongTag(vit->second, vIDTag, &(vit->first));
  }

  // Reserve tags used for model faces
  for (int i = 0; i < nBFaces ; ++i) {
    // TODO: How to assert, when bFaces is already all integers?
    minAvbl = (minAvbl <= bFaces[i][0]) ? (bFaces[i][0]+1) : minAvbl;
    PCU_ALWAYS_ASSERT(minAvbl < std::numeric_limits<int>::max());
  }

  // Set classifn tags for vertices
  for (size_t i = 0; i < mesh->count(0); ++i) {
    if (isModelVert[i]) {
      tagData[0] = 0;
      tagData[1] = minAvbl;
      minAvbl++;
      mesh->setIntTag(globalToVert[i], classifnTag, tagData);
    }
  }

  // Classification on boundary faces and their closure
  apf::Downward facesAdjToRegion, edgesAdjToFace, vertsAdjToEdge;
  ModelEdgeTags modelEdgeTags;
  int faceIdInReg = 12;
  for (int i = 0; i < nBFaces; ++i) {
    mesh->getDownward(globalToRegion[bFaces[i][1]], 2, facesAdjToRegion);
    faceIdInReg = getFaceIdInRegion(mesh,
				    globalToRegion[bFaces[i][1]], bFaces[i]);
    apf::MeshEntity* face = facesAdjToRegion[faceIdInReg];
    tagData[0] = 2;
    tagData[1] = bFaces[i][0];
    mesh->setIntTag(face, classifnTag, tagData);

    // Tag closure for classification
    mesh->getDownward(face, 1, edgesAdjToFace);
    for (int j = 0; j < 3; ++j) {
      if (mesh->hasTag(edgesAdjToFace[j], classifnTag)) {
        mesh->getIntTag(edgesAdjToFace[j], classifnTag, tagData);
        PCU_ALWAYS_ASSERT(tagData[0] == 2);
        if (tagData[0] != 2 || tagData[1] != bFaces[i][0]) {
          // Mesh edge is classified previously in another boundary
          // face's loop on a separate model face,
          // thus it has to be on a model edge
          newTagData[0] = 1;
          newTagData[1] = getModelEdgeTag(tagData[1], bFaces[i][0],
					  modelEdgeTags, &minAvbl);
          mesh->setIntTag(edgesAdjToFace[j], classifnTag, newTagData);
        } else {
          newTagData[0] = tagData[0];
          newTagData[1] = tagData[1];
        }
      } else {
        // Edge is seen first time, apply face's classification
        newTagData[0] = 2;
        newTagData[1] = bFaces[i][0];
        mesh->setIntTag(edgesAdjToFace[j], classifnTag, newTagData);
      }
      // vertices of the edge
      mesh->getDownward(edgesAdjToFace[j], 0, vertsAdjToEdge);
      for (int k = 0; k < 2; ++k) {
        if (mesh->hasTag(vertsAdjToEdge[k], classifnTag)) {
          mesh->getIntTag(vertsAdjToEdge[k], classifnTag, tagData);
          if (tagData[0] > newTagData[0]) {
            mesh->setIntTag(vertsAdjToEdge[k], classifnTag, newTagData);
          }
        } else {
          // Vertex has no classification, use the edge's
          mesh->setIntTag(vertsAdjToEdge[k], classifnTag, newTagData);
        }
      }
    }
  }

  MeshMDS* m = static_cast<MeshMDS*>(mesh);
  if ((classifnTag)) {
    int tagData[2];
    MeshEntity* ent;
    mds_id id;
    for (int dim = m->getDimension(); dim >= 0; --dim) {
      apf::MeshIterator* it = m->begin(dim);
      while ((ent = m->iterate(it))) {
        id = fromEnt(ent);
        if (m->hasTag(ent, classifnTag)) {
          m->getIntTag(ent, classifnTag, tagData);
          mds_update_model_for_entity(m->mesh, id, tagData[0], tagData[1]);
        } else {
          mds_update_model_for_entity(m->mesh, id, m->getDimension(), 
              DEFAULT_REGION_ID);
        }
      }
    }
  }

  mesh->destroyTag(classifnTag);
  mesh->destroyTag(vIDTag);
}

void derive2DMdlFromManifold(Mesh2* mesh, bool* isModelVert,
			     int nBEdges, int (*bEdges)[4],
			     GlobalToVert &globalToVert,
			     std::map<int, apf::MeshEntity*> &globalToFace)
{
  PCU_ALWAYS_ASSERT_VERBOSE(!mesh->findTag("_classifn_data"),
          "MeshTag name \"_classifn_data\" is used internally in this method\n");
  apf::MeshTag* classifnTag = mesh->createIntTag("_classifn_data", 2);
  int tagData[2], newTagData[2];
  long minAvbl = 1;
  // This is set by during apf::construct, and using anything else leads to an
  // additional region
  long DEFAULT_REGION_ID = 0;

  PCU_ALWAYS_ASSERT_VERBOSE(!mesh->findTag("_vert_id"),
          "MeshTag name \"_vert_id\" is used internally in this method\n");
  apf::MeshTag* vIDTag = mesh->createLongTag("_vert_id", 1);
  for (apf::GlobalToVert::iterator vit = globalToVert.begin();
       vit != globalToVert.end(); vit++) {
    mesh->setLongTag(vit->second, vIDTag, &(vit->first));
  }

  // Reserve tags used for model edges
  for (int i = 0; i < nBEdges; ++i) {
    // TODO: How to assert, when bEdges is already all integers?
    minAvbl = (minAvbl <= bEdges[i][0]) ? (bEdges[i][0]+1) : minAvbl;
    PCU_ALWAYS_ASSERT(minAvbl < std::numeric_limits<int>::max());
  }

  // Set classifn tags for vertices
  for (size_t i = 0; i < mesh->count(0); ++i) {
    if (isModelVert[i]) {
      tagData[0] = 0;
      tagData[1] = minAvbl;
      minAvbl++;
      mesh->setIntTag(globalToVert[i], classifnTag, tagData);
    }
  }

  // Classification of boundary edges and their closure
  apf::Downward edgesAdjToFace, vertsAdjToEdge;
  int edgeIdInFace = 12;
  for (int i = 0; i < nBEdges; ++i) {
    mesh->getDownward(globalToFace[bEdges[i][1]], 1, edgesAdjToFace);
    edgeIdInFace = getEdgeIdInFace(mesh, globalToFace[bEdges[i][1]], bEdges[i]);
    apf::MeshEntity* edge = edgesAdjToFace[edgeIdInFace];
    tagData[0] = newTagData[0] = 1;
    tagData[1] = newTagData[1] = bEdges[i][0];
    mesh->setIntTag(edge, classifnTag, tagData);

    mesh->getDownward(edge, 0, vertsAdjToEdge);
    for (int k = 0; k <2; ++k) {
      if (mesh->hasTag(vertsAdjToEdge[k], classifnTag)) {
	mesh->getIntTag(vertsAdjToEdge[k], classifnTag, tagData);
        if (tagData[0] > newTagData[0]) {
	  mesh->setIntTag(vertsAdjToEdge[k], classifnTag, newTagData);
        }
      } else {
        // Vertex has no classification, use the edge's
        mesh->setIntTag(vertsAdjToEdge[k], classifnTag, newTagData);
      }
    }
  }

  // TODO: Use classifnTag to classify
  MeshMDS* m = static_cast<MeshMDS*>(mesh);
  if ((classifnTag)) {
    int tagData[2];
    MeshEntity* ent;
    mds_id id;
    for (int dim = m->getDimension(); dim >= 0; --dim) {
      apf::MeshIterator* it = m->begin(dim);
      while ((ent = m->iterate(it))) {
        id = fromEnt(ent);
        if (m->hasTag(ent, classifnTag)) {
          m->getIntTag(ent, classifnTag, tagData);
          mds_update_model_for_entity(m->mesh, id, tagData[0], tagData[1]);
        } else {
          mds_update_model_for_entity(m->mesh, id, m->getDimension(),
              DEFAULT_REGION_ID);
        }
      }
    }
  }

  mesh->destroyTag(classifnTag);
  mesh->destroyTag(vIDTag);
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

void disownMdsModel(Mesh2* in)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  m->ownsModel = false;
}

void setMdsMatching(Mesh2* in, bool has)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  m->isMatched = has;
}

void hackMdsAdjacency(Mesh2* in, MeshEntity* up, int i, MeshEntity* down)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  mds_hack_adjacent(&m->mesh->mds, fromEnt(up), i, fromEnt(down));
}

Mesh2* loadMdsPart(gmi_model* model, const char* meshfile, pcu::PCU *PCUObj)
{
  MeshMDS* m = new MeshMDS();
  m->init(apf::getLagrange(1), PCUObj);
  m->mesh = mds_read_smb(m->getPCU()->GetCHandle(), model, meshfile, 1, m);
  m->isMatched = false;
  m->ownsModel = true;
  initResidence(m, m->getDimension());
  return m;
}

void writeMdsPart(Mesh2* in, const char* meshfile)
{
  MeshMDS* m = static_cast<MeshMDS*>(in);
  m->mesh = mds_write_smb(m->getPCU()->GetCHandle(), m->mesh, meshfile, 1, m);
}


}

extern "C" {

void mds_write_smb_meta(struct pcu_file* file, void* mesh_cpp) {
  apf::MeshMDS* m = static_cast<apf::MeshMDS*>(mesh_cpp);
  apf::save_meta(file, m);
}

void mds_read_smb_meta(struct pcu_file* file, struct mds_apf* mesh,
                       void* mesh_cpp) {
  apf::MeshMDS* m = static_cast<apf::MeshMDS*>(mesh_cpp);
/* hack warning: in order for apf::restore_data to work,
   the mds_apf pointer needs to be connected to the MeshMDS class,
   but that is typically done right after calling mds_read_smb()
   and this code is executing as a callback during read_smb() */
  m->mesh = mesh;
  apf::restore_meta(file, m);
}

}
