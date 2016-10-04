#include "apfSIM.h"
#include <apf.h>
#include <apfShape.h>
#include <SimModel.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf_simConfig.h>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#ifdef USE_FIELDSIM

#include "apfSIMDataOf.h"

apf::Field* apf::createSIMField(Mesh* m, const char* name, int valueType,
    FieldShape* shape)
{
  return makeField(m, name, valueType, 0,shape, new SIMDataOf<double>);
}

::Field* apf::getSIMField(apf::Field* f)
{
  apf::SIMDataOf<double>* data = dynamic_cast<apf::SIMDataOf<double>*>(f->getData());
  return data->getSimField();
}

apf::Field* apf::wrapSIMField(Mesh* m, pField fd)
{
  pPolyField pf = static_cast<pPolyField>(Field_def(fd));
  int order = PolyField_entOrder(pf, 0);
  apf::FieldShape* shape = apf::getLagrange(order);
  char const* name = Field_name(fd);
  int num_comp = Field_numComp(fd);
  int valueType = -1;
  switch (num_comp) {
    case 1: valueType = apf::SCALAR; break;
    case 3: valueType = apf::VECTOR; break;
    case 9: valueType = apf::MATRIX; break;
  }
  return makeField(m, name, valueType, 0, shape, new SIMDataOf<double>(fd));
}

#else

apf::Field* apf::createSIMField(Mesh* m, const char* name, int valueType,
    FieldShape* shape)
{
  (void)m;
  (void)name;
  (void)valueType;
  (void)shape;
  apf::fail("SimField not found when APF_SIM compiled");
}

::Field* apf::getSIMField(apf::Field* f)
{
  (void)f;
  apf::fail("SimField not found when APF_SIM compiled");
}

apf::Field* apf::wrapSIMField(Mesh* m, pField fd)
{
  (void)m;
  (void)fd;
  apf::fail("SimField not found when APF_SIM compiled");
}

#endif

namespace apf {

apf::Field* createSIMLagrangeField(Mesh* m, const char* name, int valueType, int order)
{
  return createSIMField(m, name, valueType, getLagrange(order));
}

apf::Field* createSIMFieldOn(Mesh* m, const char* name, int valueType)
{
  return createSIMField(m, name, valueType, m->getShape());
}

MeshSIM::MeshSIM(pParMesh m):
  mesh(m)
{
  part = PM_mesh(mesh,0);
  d = M_numRegions(part) ? 3 : 2;
  iterDim = -1;
  model = gmi_import_sim(M_model(part));
}

MeshSIM::~MeshSIM()
{
  gmi_destroy(model);
}

int MeshSIM::getDimension()
{
  return d;
}

class IteratorSIM
{
  public:
    virtual ~IteratorSIM() {}
    virtual MeshEntity* iterate() = 0;
};

class VertexIteratorSIM : public IteratorSIM
{
  public:
    VertexIteratorSIM(pMesh part)
    {iterator = M_vertexIter(part);}
    virtual ~VertexIteratorSIM()
    {VIter_delete(iterator);}
    virtual MeshEntity* iterate()
    {return reinterpret_cast<MeshEntity*>(VIter_next(iterator));}
  private:
    VIter iterator;
};

class EdgeIteratorSIM : public IteratorSIM
{
  public:
    EdgeIteratorSIM(pMesh part)
    {iterator = M_edgeIter(part);}
    virtual ~EdgeIteratorSIM()
    {EIter_delete(iterator);}
    virtual MeshEntity* iterate()
    {return reinterpret_cast<MeshEntity*>(EIter_next(iterator));}
  private:
    EIter iterator;
};

class FaceIteratorSIM : public IteratorSIM
{
  public:
    FaceIteratorSIM(pMesh part)
    {iterator = M_faceIter(part);}
    virtual ~FaceIteratorSIM()
    {FIter_delete(iterator);}
    virtual MeshEntity* iterate()
    {return reinterpret_cast<MeshEntity*>(FIter_next(iterator));}
  private:
    FIter iterator;
};

class RegionIteratorSIM : public IteratorSIM
{
  public:
    RegionIteratorSIM(pMesh part)
    {iterator = M_regionIter(part);}
    virtual ~RegionIteratorSIM()
    {RIter_delete(iterator);}
    virtual MeshEntity* iterate()
    {return reinterpret_cast<MeshEntity*>(RIter_next(iterator));}
  private:
    RIter iterator;
};

std::size_t MeshSIM::count(int dimension)
{
  if (dimension == 0)
    return M_numVertices(part);
  if (dimension == 1)
    return M_numEdges(part);
  if (dimension == 2)
    return M_numFaces(part);
  return M_numRegions(part);
}

MeshIterator* MeshSIM::begin(int dimension)
{
  IteratorSIM* iterator = 0;
  if (dimension == 0)
    iterator = new VertexIteratorSIM(part);
  if (dimension == 1)
    iterator = new EdgeIteratorSIM(part);
  if (dimension == 2)
    iterator = new FaceIteratorSIM(part);
  if (dimension == 3)
    iterator = new RegionIteratorSIM(part);
  return reinterpret_cast<MeshIterator*>(iterator);
}

MeshEntity* MeshSIM::iterate(MeshIterator* it)
{
  return reinterpret_cast<IteratorSIM*>(it)->iterate();
}

void MeshSIM::end(MeshIterator* it)
{
  delete reinterpret_cast<IteratorSIM*>(it);
}

bool MeshSIM::isShared(MeshEntity* e)
{
  return EN_isOnPartBdry(reinterpret_cast<pEntity>(e));
}

bool MeshSIM::isOwned(MeshEntity* e)
{
  return EN_isOwnerProc(reinterpret_cast<pEntity>(e));
}

int MeshSIM::getOwner(MeshEntity* e)
{
  return EN_ownerProc(reinterpret_cast<pEntity>(e));
}

static int pListToArray(pPList list, MeshEntity** array)
{
  int n = PList_size(list);
  for (int i=0; i < n; ++i)
    array[i] = reinterpret_cast<MeshEntity*>(PList_item(list,i));
  PList_delete(list);
  return n;
}

void MeshSIM::getAdjacent(MeshEntity* e,
    int dimension,
    DynamicArray<MeshEntity*>& adjacent)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  eType ent_type = EN_type(entity);
  if (apf::getDimension(this, e) == dimension)
  {
    adjacent.setSize(1);
    adjacent[0] = e;
    return;
  }
  if (ent_type == Tvertex)
  {
    pVertex vertex = static_cast<pVertex>(entity);
    if (dimension == 1)
    {
      int n = V_numEdges(vertex);
      adjacent.setSize(n);
      for (int i=0; i < n; ++i)
        adjacent[i] = reinterpret_cast<MeshEntity*>(V_edge(vertex,i));
    }
    if (dimension == 2)
      pListToDynamicArray<MeshEntity*>(V_faces(vertex),adjacent);
    if (dimension == 3)
      pListToDynamicArray<MeshEntity*>(V_regions(vertex),adjacent);
  }
  if (ent_type == Tedge)
  {
    pEdge edge = static_cast<pEdge>(entity);
    if (dimension == 0)
    {
      adjacent.setSize(2);
      for (int i=0; i < 2; ++i)
        adjacent[i] = reinterpret_cast<MeshEntity*>(E_vertex(edge,i));
    }
    if (dimension == 2)
    {
      int n = E_numFaces(edge);
      adjacent.setSize(n);
      for (int i=0; i < n; ++i)
        adjacent[i] = reinterpret_cast<MeshEntity*>(E_face(edge,i));
    }
    if (dimension == 3)
      pListToDynamicArray<MeshEntity*>(E_regions(edge),adjacent);
  }
  if (ent_type == Tface)
  {
    pFace face = static_cast<pFace>(entity);
    if (dimension == 0)
      pListToDynamicArray<MeshEntity*>(F_vertices(face,1),adjacent);
    if (dimension == 1)
    {
      int n = F_numEdges(face);
      adjacent.setSize(n);
      for (int i=0; i < n; ++i)
        adjacent[i] = reinterpret_cast<MeshEntity*>(F_edge(face,i));
    }
    if (dimension == 3)
    {
      adjacent.setSize(countUpward(e));
      int a=0;
      for (int i=0; i < 2; ++i)
      {
        MeshEntity * me =  reinterpret_cast<MeshEntity*>(F_region(face,i));
        if(me != NULL)
          adjacent[a++] = me;
      }
    }
  }
  if (ent_type == Tregion)
  {
    pRegion region = static_cast<pRegion>(entity);
    if (dimension == 0)
      pListToDynamicArray<MeshEntity*>(R_vertices(region,1),adjacent);
    if (dimension == 1)
      pListToDynamicArray<MeshEntity*>(R_edges(region,1),adjacent);
    if (dimension == 2)
    {
      int n = R_numFaces(region);
      adjacent.setSize(n);
      for (int i=0; i < n; ++i)
        adjacent[i] = reinterpret_cast<MeshEntity*>(R_face(region,i));
    }
  }
}

int MeshSIM::getDownward(MeshEntity* e,
    int dimension,
    MeshEntity** down)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  eType ent_type = EN_type(entity);
  if (apf::getDimension(this, e) == dimension)
  {
    down[0] = e;
  }
  else if (ent_type == Tedge)
  {
    pEdge edge = static_cast<pEdge>(entity);
    for (int i=0; i < 2; ++i)
      down[i] = reinterpret_cast<MeshEntity*>(E_vertex(edge,i));
  }
  else if (ent_type == Tface)
  {
    pFace face = static_cast<pFace>(entity);
    if (dimension == 0)
      pListToArray(F_vertices(face,1),down);
    if (dimension == 1)
    {
      int n = F_numEdges(face);
      for (int i=0; i < n; ++i)
        down[i] = reinterpret_cast<MeshEntity*>(F_edge(face,i));
    }
  }
  else if (ent_type == Tregion)
  {
    pRegion region = static_cast<pRegion>(entity);
    if (dimension == 0)
      pListToArray(R_vertices(region,1),down);
    if (dimension == 1)
      pListToArray(R_edges(region,1),down);
    if (dimension == 2)
    {
      int n = R_numFaces(region);
      for (int i=0; i < n; ++i)
        down[i] = reinterpret_cast<MeshEntity*>(R_face(region,i));
    }
  }
  return Mesh::adjacentCount[getType(e)][dimension];
}

int MeshSIM::countUpward(MeshEntity* e)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  eType ent_type = EN_type(entity);
  if (ent_type == Tvertex)
    return V_numEdges(reinterpret_cast<pVertex>(entity));
  if (ent_type == Tedge)
    return E_numFaces(reinterpret_cast<pEdge>(entity));
  if (ent_type == Tface)
  {
    int n = 0;
    for (int i=0; i < 2; ++i)
      if (F_region(reinterpret_cast<pFace>(entity),i))
        ++n;
    return n;
  }
  return 0;
}

MeshEntity* MeshSIM::getUpward(MeshEntity* e, int i)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  eType ent_type = EN_type(entity);
  if (ent_type == Tvertex)
    return reinterpret_cast<MeshEntity*>(
        V_edge(reinterpret_cast<pVertex>(entity),i));
  if (ent_type == Tedge)
    return reinterpret_cast<MeshEntity*>(
        E_face(reinterpret_cast<pEdge>(entity),i));
  if (ent_type == Tface)
  {
    int a = 0;
    for (int i=0; i < 2; ++i)
    {
      MeshEntity* r = reinterpret_cast<MeshEntity*>(
          F_region(reinterpret_cast<pFace>(entity),i));
      if ((r)&&(a++ == i)) return r;
    }
  }
  return 0;
}

void MeshSIM::getUp(MeshEntity* e, Up& up)
{
  up.n = countUpward(e);
  for (int i = 0; i < up.n; ++i)
    up.e[i] = getUpward(e,i);
}

bool MeshSIM::hasUp(MeshEntity* e)
{
  return countUpward(e) != 0;
}

void MeshSIM::getPoint_(MeshEntity* e, int node, Vector3& point)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  eType type = EN_type(entity);
  pPoint pt = 0;
  if (type == Tvertex)
    pt = V_point(static_cast<pVertex>(entity));
  if (type == Tedge)
    pt = E_point(static_cast<pEdge>(entity),node);
  P_coord(pt,&(point[0]));
}

void MeshSIM::getParam(MeshEntity* e, Vector3& point)
{
  pVertex v = reinterpret_cast<pVertex>(e);
  pPoint pt = V_point(v);
  int d = getModelType(toModel(e));
  if (d==1)
    point[0] = P_param1(pt);
  if (d==2)
  {
    int obsolete_patch;
    P_param2(pt,&(point[0]),&(point[1]),&obsolete_patch);
  }
}

Mesh::Type MeshSIM::getType(MeshEntity* e)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  eType ent_type = EN_type(entity);
  if(ent_type == Tvertex)
    return Mesh::VERTEX;
  else if(ent_type == Tedge)
    return Mesh::EDGE;
  else if(ent_type == Tface)
  {
    int num_edges = F_numEdges(static_cast<pFace>(entity));
    if(num_edges == 3)
      return Mesh::TRIANGLE;
    else if(num_edges == 4)
      return Mesh::QUAD;
  }
  else
  {
    rType reg_type = R_topoType(static_cast<pRegion>(entity));
    if(reg_type == Rtet)
      return Mesh::TET;
    else if(reg_type == Rwedge)
      return Mesh::PRISM;
    else if(reg_type == Rpyramid)
      return Mesh::PYRAMID;
    else if(reg_type == Rhex)
      return Mesh::HEX;
  }
  abort();
}

void MeshSIM::getRemotes(MeshEntity* e, Copies& remotes)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  pEntCopies copies = EN_copies(entity);
  if (!copies) return;
  for (int i=0; i < EntCopies_size(copies); ++i)
    remotes[PMU_proc(EntCopies_gid(copies,i))] = reinterpret_cast<MeshEntity*>(
        EntCopies_ent(copies,i));
}

void MeshSIM::getResidence(MeshEntity* e, Parts& residence)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  pEntCopies copies = EN_copies(entity);
  residence.insert(getId());
  if (!copies) return;
  for (int i=0; i < EntCopies_size(copies); ++i)
    residence.insert(PMU_proc(EntCopies_gid(copies,i)));
}

class TagSIM
{
  public:
    TagSIM(pParMesh m,
           const char* n,
           std::size_t unitSize,
           int c):
      count(c),
      mesh(m),
      name(n)
    {
      id = MD_newMeshDataId(n);
      comm = PM_newAttachDataCommu(unitSize/sizeof(int),0,count);
      PM_setMigrId(mesh,id);
    }
    virtual ~TagSIM()
    {
      MD_removeMeshCallback(id,CBdelete);
      MD_removeMeshCallback(id,CBmigrateOut);
      MD_removeMeshCallback(id,CBmigrateIn);
      PM_removeMigrId(mesh,id);
      MD_deleteMeshDataId(id);
      AttachDataCommu_delete(comm);
    }
    virtual void* allocate() = 0;
    virtual void deallocate(void* p) = 0;
    virtual int getType() = 0;
    bool has(MeshEntity* e)
    {
      pEntity entity = reinterpret_cast<pEntity>(e);
      return EN_getDataPtr(entity,id,NULL);
    }
    void set(MeshEntity* e, void* p)
    {
      pEntity entity = reinterpret_cast<pEntity>(e);
      EN_attachDataPtr(entity,id,p);
    }
    void* get(MeshEntity* e)
    {
      pEntity entity = reinterpret_cast<pEntity>(e);
      if ( ! has(e))
        set(e,this->allocate());
      void* p;
      EN_getDataPtr(entity,id,&p);
      return p;
    }
    void remove(MeshEntity* e)
    {
      pEntity entity = reinterpret_cast<pEntity>(e);
      this->deallocate(this->get(e));
      EN_deleteData(entity,id);
    }
    int count;
    pParMesh mesh;
    pMeshDataId id;
    pAttachDataCommu comm;
    std::string name; //Simmetrix has no "get tag name" API
};

class DoubleTagSIM : public TagSIM
{
  public:
    static int deleteDoubleCB(
        void*,
        pAttachDataId id,
        int ev,
        void** data,
        void*)
    {
      static_cast<DoubleTagSIM*>(MD_callbackData(static_cast<pMeshDataId>(id), ev))->deallocate(*data);
      *data = 0;
      return 1;
    }

    DoubleTagSIM(pParMesh m, const char* name, int c):
      TagSIM(m,name,sizeof(double),c)
    {
      MD_setMeshCallback(id,CBdelete,deleteDoubleCB,this);
      MD_setMeshCallback(id,CBmigrateOut,pm_sendDblArray,comm);
      MD_setMeshCallback(id,CBmigrateIn,pm_recvDblArray,comm);
    }
    virtual void* allocate()
    {
      return count == 1 ? new double() : new double[count]();
    }
    virtual void deallocate(void* p)
    {
      double* p2 = static_cast<double*>(p);
      if (count == 1)
        delete p2;
      else
        delete [] p2;
    }
    virtual int getType() {return Mesh::DOUBLE;}
    void get(MeshEntity* e, double* p)
    {
      double* internal = static_cast<double*>(TagSIM::get(e));
      for (int i=0; i < count; ++i)
        p[i] = internal[i];
    }
    void set(MeshEntity* e, double const* p)
    {
      double* internal = static_cast<double*>(TagSIM::get(e));
      for (int i=0; i < count; ++i)
        internal[i] = p[i];
    }
};

class IntTagSIM : public TagSIM
{
  public:
    static int deleteIntCB(
        void*,
        pAttachDataId id,
        int ev,
        void** data,
        void*)
    {
      static_cast<IntTagSIM*>(MD_callbackData(static_cast<pMeshDataId>(id), ev))->deallocate(*data);
      *data = 0;
      return 1;
    }

    IntTagSIM(pParMesh m, const char* name, int c):
      TagSIM(m,name,sizeof(int),c)
    {
      MD_setMeshCallback(id,CBdelete,deleteIntCB,this);
      MD_setMeshCallback(id,CBmigrateOut,pm_sendIntArray,comm);
      MD_setMeshCallback(id,CBmigrateIn,pm_recvIntArray,comm);
    }
    virtual void* allocate()
    {
      return count == 1 ? new int() : new int[count]();
    }
    virtual void deallocate(void* p)
    {
      int* p2 = static_cast<int*>(p);
      if (count == 1)
        delete p2;
      else
        delete [] p2;
    }
    virtual int getType() {return Mesh::INT;}
    void get(MeshEntity* e, int* p)
    {
      int* internal = static_cast<int*>(TagSIM::get(e));
      for (int i=0; i < count; ++i)
        p[i] = internal[i];
    }
    void set(MeshEntity* e, int const* p)
    {
      int* internal = static_cast<int*>(TagSIM::get(e));
      for (int i=0; i < count; ++i)
        internal[i] = p[i];
    }
};

class LongTagSIM : public TagSIM
{
  public:
    static int deleteLongCB(
        void*,
        pAttachDataId id,
        int ev,
        void** data,
        void*)
    {
      static_cast<LongTagSIM*>(MD_callbackData(static_cast<pMeshDataId>(id), ev))->deallocate(*data);
      *data = 0;
      return 1;
    }

    LongTagSIM(pParMesh m, const char* name, int c):
      TagSIM(m,name,sizeof(long),c)
    {
      MD_setMeshCallback(id,CBdelete,deleteLongCB,this);
      /* note: long tags won't get auto-migration support until
         this is filled in: */
      // MD_setMeshCallback(id,CBmigrateOut,pm_sendIntArray,comm);
      // MD_setMeshCallback(id,CBmigrateIn,pm_recvIntArray,comm);
      PM_removeMigrId(m, id);
    }
    virtual void* allocate()
    {
      return count == 1 ? new long() : new long[count]();
    }
    virtual void deallocate(void* p)
    {
      long* p2 = static_cast<long*>(p);
      if (count == 1)
        delete p2;
      else
        delete [] p2;
    }
    virtual int getType() {return Mesh::LONG;}
    void get(MeshEntity* e, long* p)
    {
      long* internal = static_cast<long*>(TagSIM::get(e));
      for (int i=0; i < count; ++i)
        p[i] = internal[i];
    }
    void set(MeshEntity* e, long const* p)
    {
      long* internal = static_cast<long*>(TagSIM::get(e));
      for (int i=0; i < count; ++i)
        internal[i] = p[i];
    }
};

MeshTag* MeshSIM::createDoubleTag(const char* name, int size)
{
  TagSIM* tag = new DoubleTagSIM(mesh,name,size);
  tags.push_back(tag);
  return reinterpret_cast<MeshTag*>(tag);
}

MeshTag* MeshSIM::createIntTag(const char* name, int size)
{
  TagSIM* tag = new IntTagSIM(mesh,name,size);
  tags.push_back(tag);
  return reinterpret_cast<MeshTag*>(tag);
}

MeshTag* MeshSIM::createLongTag(const char* name, int size)
{
  TagSIM* tag = new LongTagSIM(mesh,name,size);
  tags.push_back(tag);
  return reinterpret_cast<MeshTag*>(tag);
}

MeshTag* MeshSIM::findTag(const char* name)
{
  for (size_t i=0; i < tags.size(); ++i)
    if (tags[i]->name == name)
      return reinterpret_cast<MeshTag*>(tags[i]);
  return 0;
}

void MeshSIM::destroyTag(MeshTag* tag)
{
  TagSIM* tagSim = reinterpret_cast<TagSIM*>(tag);
  tags.erase(std::find(tags.begin(),tags.end(),tagSim));
  delete tagSim;
}

void MeshSIM::renameTag(MeshTag*, const char*)
{
  apf::fail("MeshSIM::renameTag called!\n");
}

void MeshSIM::getTags(DynamicArray<MeshTag*>& ts)
{
  ts.setSize(tags.size());
  for (size_t i=0; i < tags.size(); ++i)
    ts[i] = reinterpret_cast<MeshTag*>(tags[i]);
}

void MeshSIM::getDoubleTag(MeshEntity* e, MeshTag* tag, double* data)
{
  DoubleTagSIM* tagSim = reinterpret_cast<DoubleTagSIM*>(tag);
  tagSim->get(e,data);
}

void MeshSIM::setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data)
{
  DoubleTagSIM* tagSim = reinterpret_cast<DoubleTagSIM*>(tag);
  tagSim->set(e,data);
}

void MeshSIM::getIntTag(MeshEntity* e, MeshTag* tag, int* data)
{
  IntTagSIM* tagSim = reinterpret_cast<IntTagSIM*>(tag);
  tagSim->get(e,data);
}

void MeshSIM::setIntTag(MeshEntity* e, MeshTag* tag, int const* data)
{
  IntTagSIM* tagSim = reinterpret_cast<IntTagSIM*>(tag);
  tagSim->set(e,data);
}

void MeshSIM::getLongTag(MeshEntity* e, MeshTag* tag, long* data)
{
  LongTagSIM* tagSim = reinterpret_cast<LongTagSIM*>(tag);
  tagSim->get(e,data);
}

void MeshSIM::setLongTag(MeshEntity* e, MeshTag* tag, long const* data)
{
  LongTagSIM* tagSim = reinterpret_cast<LongTagSIM*>(tag);
  tagSim->set(e,data);
}

void MeshSIM::removeTag(MeshEntity* e, MeshTag* tag)
{
  TagSIM* tagSim = reinterpret_cast<TagSIM*>(tag);
  tagSim->remove(e);
}

bool MeshSIM::hasTag(MeshEntity* e, MeshTag* tag)
{
  TagSIM* tagSim = reinterpret_cast<TagSIM*>(tag);
  return tagSim->has(e);
}

int MeshSIM::getTagType(MeshTag* t)
{
  TagSIM* tag = reinterpret_cast<TagSIM*>(t);
  return tag->getType();
}

int MeshSIM::getTagSize(MeshTag* t)
{
  TagSIM* tag = reinterpret_cast<TagSIM*>(t);
  return tag->count;
}

const char* MeshSIM::getTagName(MeshTag* t)
{
  TagSIM* tag = reinterpret_cast<TagSIM*>(t);
  return tag->name.c_str();
}

ModelEntity* MeshSIM::toModel(MeshEntity* e)
{
  pEntity entity = reinterpret_cast<pEntity>(e);
  return reinterpret_cast<ModelEntity*>(EN_whatIn(entity));
}

gmi_model* MeshSIM::getModel()
{
  return model;
}

void MeshSIM::migrate(Migration* plan)
{
  pMigrator migrator = PM_newMigrator(mesh,sthreadNone);
  Migrator_reset(migrator,this->getDimension());
  int gid = PMU_gid(getId(),0);
  for (int i=0; i < plan->count(); ++i)
  {
    MeshEntity* e = plan->get(i);
    int newgid = PMU_gid(plan->sending(e),0);
    pEntity entity = reinterpret_cast<pEntity>(e);
    Migrator_add(migrator,entity,newgid,gid);
  }
  delete plan;
  Migrator_run(migrator,NULL);
  Migrator_delete(migrator);
}

int MeshSIM::getId()
{
  return PMU_gid(PMU_rank(),0);
}

void MeshSIM::writeNative(const char* fileName)
{
  PM_write(mesh,fileName,sthreadNone,NULL);
}

void MeshSIM::destroyNative()
{
  M_release(mesh);
}

void MeshSIM::verify()
{
  assert(PM_verify(mesh,0,sthreadNone,NULL));
}

void MeshSIM::getMatches(MeshEntity* e, Matches& m)
{
  pEntity ent = reinterpret_cast<pEntity>(e);
  pPList l = EN_getMatchingEnts(ent, NULL);
  if (!l)
    return;
  int n = PList_size(l);
  m.setSize(n - 1);
  int j = 0;
  for (int i=0; i < n; ++i)
  {
    MeshEntity* match = reinterpret_cast<MeshEntity*>(PList_item(l,i));
    if (match == e)
      continue;
    m[j].peer = this->getId();
    m[j].entity = match;
    j++;
  }
  assert(j == n - 1);
  PList_delete(l);
}

void MeshSIM::setPoint_(MeshEntity * me, int node, Vector3 const & p)
{
  pEntity entity = reinterpret_cast<pEntity>(me);
  eType type = EN_type(entity);
  pPoint pt = 0;
  if (type == Tvertex)
    pt = V_point(static_cast<pVertex>(entity));
  if (type == Tedge)
    pt = E_point(static_cast<pEdge>(entity),node);
  P_setPos(pt,p[0],p[1],p[2]);
}

static bool isQuadratic(pParMesh mesh)
{
  pMesh part = PM_mesh(mesh,0);
  EIter it = M_edgeIter(part);
  pEdge e;
  bool result = true;
  while ((e = EIter_next(it)))
  {
    if (E_numPoints(e) != 1)
    {
      result = false;
      break;
    }
  }
  EIter_delete(it);
  return result;
}

static bool findMatches(Mesh* m)
{
  bool found = false;
  for (int i = 0; i < 4; ++i)
  {
    MeshIterator* it = m->begin(i);
    MeshEntity* e;
    while ((e = m->iterate(it)))
    {
      pEntity ent = reinterpret_cast<pEntity>(e);
      pPList l = EN_getMatchingEnts(ent, NULL);
      if (l)
      {
        found = true;
        PList_delete(l);
        break;
      }
    }
    m->end(it);
    if (found)
      break;
  }
  return found;
}

Mesh2* createMesh(pParMesh mesh)
{
  /* require one part per process currently for SIM */
  assert(PM_numParts(mesh)==1);
  MeshSIM* m = new MeshSIM(mesh);
  int order = 1;
  if (isQuadratic(mesh))
    order = 2;
  m->init(getLagrange(order));
  m->hasMatches = findMatches(m);
  return m;
}

MeshEntity* castEntity(pEntity entity)
{
  return reinterpret_cast<MeshEntity*>(entity);
}

}//namespace apf
