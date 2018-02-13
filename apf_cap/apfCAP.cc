#include "apfCAP.h"
#include <PCU.h>
#include <apf.h>
#include <apfShape.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <cstdlib>
#include <pcu_util.h>
#include <algorithm>


namespace apf {


MeshEntity* toEntity(M_MTopo topo)
{
  std::size_t hdl = topo.get();
  PCU_ALWAYS_ASSERT(hdl > 0);
  return reinterpret_cast<MeshEntity*>(hdl);
}

M_MTopo fromEntity(MeshEntity* e)
{
  std::size_t hdl = reinterpret_cast<std::size_t>(e);
  M_MTopo topo;
  topo.set(hdl);
  return topo;
}

/* static void setupAdjacencies(MeshDatabaseInterface* mdb) */
/* { */
/*   // downward adjacencies */
/*   MStatus st; */
/*   if( !mdb->adjacency_exists(TOPO_EDGE, TOPO_VERTEX)) { */
/*     st = mdb->compute_adjacency(TOPO_EDGE, TOPO_VERTEX); */
/*     PCU_ALWAYS_ASSERT(st == STATUS_OK); */
/*   } */

/*   if( !mdb->adjacency_exists(TOPO_FACE, TOPO_VERTEX)) { */
/*     st = mdb->compute_adjacency(TOPO_FACE, TOPO_VERTEX); */
/*     PCU_ALWAYS_ASSERT(st == STATUS_OK); */
/*   } */
/*   if( !mdb->adjacency_exists(TOPO_FACE, TOPO_EDGE)) { */
/*     st = mdb->compute_adjacency(TOPO_FACE, TOPO_EDGE); */
/*     PCU_ALWAYS_ASSERT(st == STATUS_OK); */
/*   } */

/*   if( !mdb->adjacency_exists(TOPO_REGION, TOPO_VERTEX)) { */
/*     st = mdb->compute_adjacency(TOPO_REGION, TOPO_VERTEX); */
/*     PCU_ALWAYS_ASSERT(st == STATUS_OK); */
/*   } */
/*   if( !mdb->adjacency_exists(TOPO_REGION, TOPO_EDGE)) { */
/*     st = mdb->compute_adjacency(TOPO_REGION, TOPO_EDGE); */
/*     PCU_ALWAYS_ASSERT(st == STATUS_OK); */
/*   } */
/*   if( !mdb->adjacency_exists(TOPO_REGION, TOPO_FACE)) { */
/*     st = mdb->compute_adjacency(TOPO_REGION, TOPO_FACE); */
/*     PCU_ALWAYS_ASSERT(st == STATUS_OK); */
/*   } */
/*   // upward adjacencies */
/* } */

MeshCAP::MeshCAP(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb):
  meshInterface(mdb), geomInterface(gdb)
{
  PCU_ALWAYS_ASSERT(meshInterface);

  std::size_t numRegions = 0;
  meshInterface->get_num_topos(TOPO_REGION, numRegions);
  d = numRegions ? 3 : 2;
  iterDim = -1;
  model = gmi_import_cap(geomInterface);
  /* setupAdjacencies(meshInterface); */
  /* MStatus stat = meshInterface->compute_adjacency(); */
  /* /1* MStatus stat = meshInterface->compute_adjacency(TOPO_FACE, TOPO_EDGE); *1/ */
  /* PCU_ALWAYS_ASSERT(stat == STATUS_OK); */
}

MeshCAP::~MeshCAP()
{
}

int MeshCAP::getDimension()
{
  return d;
}

std::size_t MeshCAP::count(int dimension)
{
  std::size_t count = 0;
  if (dimension == 0)
    meshInterface->get_num_topos(TOPO_VERTEX, count);
  if (dimension == 1)
    meshInterface->get_num_topos(TOPO_EDGE, count);
  if (dimension == 2)
    meshInterface->get_num_topos(TOPO_FACE, count);
  if (dimension == 3)
    meshInterface->get_num_topos(TOPO_REGION, count);
  return count;
}

Mesh::Type MeshCAP::getType(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  MeshShape topoShape;
  meshInterface->get_topo_shape(topo, topoShape);
  if (topoShape == SHAPE_NODE)
    return Mesh::VERTEX;
  else if (topoShape == SHAPE_SEGMENT)
    return Mesh::EDGE;
  else if (topoShape == SHAPE_TRIANGLE)
    return Mesh::TRIANGLE;
  else if (topoShape == SHAPE_QUAD)
    return Mesh::QUAD;
  else if (topoShape == SHAPE_TETRA)
    return Mesh::TET;
  else if (topoShape == SHAPE_HEX)
    return Mesh::HEX;
  else if (topoShape == SHAPE_PRISM)
    return Mesh::PRISM;
  else if (topoShape == SHAPE_PYRAMID)
    return Mesh::PYRAMID;
  else
    apf::fail("MeshCAP::getType encountered an unknown entity type!\n");
}

void MeshCAP::verify()
{
  apf::fail("MeshCAP::verify called!\n");
}

void MeshCAP::writeNative(const char* fileName)
{
  (void)fileName;
  apf::fail("MeshCAP::writeNative called!\n");
}

void MeshCAP::destroyNative()
{
  apf::fail("MeshCAP::destroyNative called!\n");
}

void MeshCAP::getPoint_(MeshEntity* e, int node, Vector3& point)
{
  (void)node;
  M_MTopo topo = fromEntity(e);
  if (meshInterface->is_vertex(topo))
    meshInterface->get_vertex_coord(topo, &(point[0]));
  else
    apf::fail("MeshCAP::getPoint_ is called for entity other than vertex!\n");
}

void MeshCAP::setPoint_(MeshEntity * me, int node, Vector3 const & p)
{
  (void)node;
  M_MTopo topo = fromEntity(me);
  if (meshInterface->is_vertex(topo))
    meshInterface->set_vertex_coord(topo, &(p[0]));
  else
    apf::fail("MeshCAP::getPoint_ is called for entity other than vertex!\n");
}

void MeshCAP::getParam(MeshEntity* e, Vector3& point)
{
  M_MTopo topo = fromEntity(e);
  /* int d = getModelType(toModel(e)); */
  /* PCU_ALWAYS_ASSERT(d==1 || d==2); */
  double u, v;
  GeometryTopoType gtype;
  meshInterface->get_vertex_uv_parameters(topo, u, v, gtype);
  point = Vector3(u, v, 0.);
}

void MeshCAP::setParam(MeshEntity* e, Vector3 const& point)
{
  M_MTopo topo = fromEntity(e);
  /* int d = getModelType(toModel(e)); */
  /* PCU_ALWAYS_ASSERT(d==1 || d==2); */
  meshInterface->set_vertex_uv_parameters(topo, point[0], point[1]);
}


MeshIterator* MeshCAP::begin(int dimension)
{
  MeshSmartIterator* miter = new MeshSmartIterator(meshInterface);
  if (dimension == 0)
    meshInterface->get_topo_iterator(TOPO_VERTEX, *miter);
  if (dimension == 1)
    meshInterface->get_topo_iterator(TOPO_EDGE, *miter);
  if (dimension == 2)
    meshInterface->get_topo_iterator(TOPO_FACE, *miter);
  if (dimension == 3)
    meshInterface->get_topo_iterator(TOPO_REGION, *miter);
  meshInterface->iterator_begin(*miter);
  return reinterpret_cast<MeshIterator*>(miter);
}

/* NOTE: miter is located at the first item in the list, therefore
 * iterate has to return it before calling iterator_next on miter
 */
MeshEntity* MeshCAP::iterate(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);

  M_MTopo topo = meshInterface->iterator_value(*miter);

  if (!meshInterface->iterator_end(*miter))
    meshInterface->iterator_next(*miter);
  else
    return 0;

  return toEntity(topo);
}

void MeshCAP::end(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);
  delete miter;
}

void MeshCAP::increment(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);
  meshInterface->iterator_next(*miter);
  it = reinterpret_cast<MeshIterator*>(miter);
}

bool MeshCAP::isDone(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);
  return meshInterface->iterator_end(*miter);
}

MeshEntity* MeshCAP::deref(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);
  M_MTopo topo = meshInterface->iterator_value(*miter);
  return toEntity(topo);
}

void MeshCAP::getAdjacent(MeshEntity* e,
    int dimension,
    DynamicArray<MeshEntity*>& adjacent)
{
  M_MTopo topo = fromEntity(e);
  MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<M_MTopo> adjTopos;
  if (apf::getDimension(this, e) == dimension)
  {
    adjacent.setSize(1);
    adjacent[0] = e;
    return;
  }
  if (type == TOPO_VERTEX)
  {
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);
    if (dimension == 3)
      meshInterface->get_adjacency_vector(topo, TOPO_REGION, adjTopos);
  }
  if (type == TOPO_EDGE)
  {
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);
    if (dimension == 3)
      meshInterface->get_adjacency_vector(topo, TOPO_REGION, adjTopos);
  }
  if (type == TOPO_FACE)
  {
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
    if (dimension == 3)
      meshInterface->get_adjacency_vector(topo, TOPO_REGION, adjTopos);
  }
  if (type == TOPO_REGION)
  {
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);
  }
  adjacent.setSize(adjTopos.size());
  for (size_t i = 0; i < adjTopos.size(); i++) {
    adjacent[i] = toEntity(adjTopos[i]);
  }
}

int MeshCAP::getDownward(MeshEntity* e,
    int dimension,
    MeshEntity** down)
{
  M_MTopo topo = fromEntity(e);
  MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<M_MTopo> adjTopos;
  if (apf::getDimension(this, e) == dimension)
  {
    down[0] = e;
    return 1;
  }
  else if (type == TOPO_EDGE)
  {
    PCU_ALWAYS_ASSERT(dimension == 0);
    meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
  }
  else if (type == TOPO_FACE)
  {
    PCU_ALWAYS_ASSERT(dimension == 0 || dimension == 1);
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
  }
  else if (type == TOPO_REGION)
  {
    PCU_ALWAYS_ASSERT(dimension == 0 || dimension == 1 || dimension == 2);
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);
  }
  for (std::size_t i = 0; i < adjTopos.size(); i++)
    down[i] = toEntity(adjTopos[i]);
  /* std::cout << adjTopos.size() << "," << Mesh::adjacentCount[getType(e)][dimension] << std::endl; */
  PCU_ALWAYS_ASSERT(adjTopos.size() == (std::size_t)Mesh::adjacentCount[getType(e)][dimension]);
  return adjTopos.size();
}

MeshEntity* MeshCAP::getUpward(MeshEntity* e, int i)
{
  M_MTopo topo = fromEntity(e);
  MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<M_MTopo> adjTopos;
  if (type == TOPO_VERTEX)
  {
    meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
    return toEntity(adjTopos[i]);
  }
  if (type == TOPO_EDGE)
  {
    meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);
    return toEntity(adjTopos[i]);
  }
  if (type == TOPO_FACE)
  {
    meshInterface->get_adjacency_vector(topo, TOPO_REGION, adjTopos);
    return toEntity(adjTopos[i]);
  }
  return 0;
}

bool MeshCAP::hasUp(MeshEntity* e)
{
  return countUpward(e) != 0;
}

bool MeshCAP::hasAdjacency(int from_dim, int to_dim)
{
  return (abs(from_dim - to_dim) == 1);
}

void MeshCAP::createAdjacency(int from_dim, int to_dim)
{
  (void)from_dim;
  (void)to_dim;
  apf::fail("MeshCAP::createAdjacency called!\n");
}

void MeshCAP::deleteAdjacency(int from_dim, int to_dim)
{
  (void)from_dim;
  (void)to_dim;
  apf::fail("MeshCAP::deleteAdjacency called!\n");
}

void MeshCAP::getUp(MeshEntity* e, Up& up)
{
  up.n = countUpward(e);
  for (int i = 0; i < up.n; i++) {
    up.e[i] = getUpward(e, i);
  }
}

int MeshCAP::countUpward(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<M_MTopo> adjTopos;
  if (type == TOPO_VERTEX)
    meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
  if (type == TOPO_EDGE)
    meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);
  if (type == TOPO_FACE)
    meshInterface->get_adjacency_vector(topo, TOPO_REGION, adjTopos);
  return (int)adjTopos.size();
}

ModelEntity* MeshCAP::toModel(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  M_GTopo gtopo;
  GeometryTopoType gtype;
  meshInterface->get_geom_entity(topo, gtype, gtopo);
  gmi_ent* g = toGmiEntity(gtopo);
  return reinterpret_cast<ModelEntity*>(g);
}

gmi_model* MeshCAP::getModel()
{
  return model;
}

void MeshCAP::setModelEntity(MeshEntity* e, ModelEntity* me)
{
  (void)e;
  (void)me;
  apf::fail("MeshCAP::setModelEntity called!\n");
}

MeshEntity* MeshCAP::createVert_(ModelEntity* me)
{
  double xyz[3] = {0., 0., 0.}; // this will be set later by setPoint_
  M_MTopo vertex; // to be created
  if ( !me ) {
    meshInterface->create_vertex(xyz, vertex);
    return toEntity(vertex);
  }
  gmi_ent* g = reinterpret_cast<gmi_ent*>(me);
  M_GTopo gtopo = fromGmiEntity(g);
  GeometryTopoType gtype;
  int d = getModelType(me);
  switch (d) {
    case 0:
      gtype = GVERTEX;
      break;
    case 1:
      gtype = GEDGE;
      break;
    case 2:
      gtype = GFACE;
      break;
    case 3:
      gtype = GREGION;
      break;
    default:
      break;
  }
  meshInterface->create_vertex(xyz, vertex, gtype, gtopo);
  return toEntity(vertex);
}

int const degree[Mesh::TYPES][4] =
{{1, 0,0,0} /* MDS_VERTEX */
,{2, 1,0,0} /* MDS_EDGE */
,{3, 3,1,0} /* MDS_TRIANGLE */
,{4, 4,1,0} /* MDS_QUADRILATERAL */
,{4, 6,4,1} /* MDS_TETRAHEDRON */
,{8,12,6,1} /* MDS_HEXAHEDRON */
,{6, 9,5,1} /* MDS_WEDGE */
,{5, 8,5,1} /* MDS_PYRAMID */
};

MeshEntity* MeshCAP::createEntity_(int type, ModelEntity* me, MeshEntity** down)
{
  MeshShape shape;
  switch (type) {
    case Mesh::VERTEX:
      shape = SHAPE_NODE;
      break;
    case Mesh::EDGE:
      shape = SHAPE_SEGMENT;
      break;
    case Mesh::TRIANGLE:
      shape = SHAPE_TRIANGLE;
      break;
    case Mesh::QUAD:
      shape = SHAPE_QUAD;
      break;
    case Mesh::TET:
      shape = SHAPE_TETRA;
      break;
    case Mesh::HEX:
      shape = SHAPE_HEX;
      break;
    case Mesh::PRISM:
      shape = SHAPE_PRISM;
      break;
    case Mesh::PYRAMID:
      shape = SHAPE_PYRAMID;
      break;
    default:
      break;
  }


  int dim = apf::getDimension(this, down[0]);
  int n = degree[type][dim];

  if (dim == 2) {
    for (int i = 0; i < 4; i++) {
      MeshEntity* dv[3];
      getDownward(down[i], 0, dv);
      for (int j = 0; j < 3; j++) {
      	Vector3 p;
        getPoint_(dv[j], 0, p);
        std::cout << p << std::endl;
      }

    }
  }

  std::vector<M_MTopo> mtopos;
  for (int i = 0; i < n; i++)
    mtopos.push_back(fromEntity(down[i]));

  M_MTopo topo;

  // first take care of the case where model is 0
  if ( !me ) {
    meshInterface->create_topo(shape, mtopos, topo);
    return toEntity(topo);
  }

  GeometryTopoType gtype;
  int d = getModelType(me);
  switch (d) {
    case 0:
      gtype = GVERTEX;
      break;
    case 1:
      gtype = GEDGE;
      break;
    case 2:
      gtype = GFACE;
      break;
    case 3:
      gtype = GREGION;
      break;
    default:
      break;
  }

  gmi_ent* g = reinterpret_cast<gmi_ent*>(me);
  M_GTopo gtopo = fromGmiEntity(g);

  meshInterface->create_topo(shape, mtopos, topo, gtype, gtopo);
  return toEntity(topo);
}

void MeshCAP::destroy_(MeshEntity* e)
{
  // remove the tags manually first
  for (std::size_t i = 0; i < tags.size(); i++)
    removeTag(e, reinterpret_cast<MeshTag*>(tags[i]));
  // now delete the entity
  M_MTopo topo = fromEntity(e);
  meshInterface->delete_topo(topo);
}

class TagCAP
{
  public:
    TagCAP(MeshDatabaseInterface* m,
	   const char* n,
           int c):
      mesh(m),
      count(c),
      name(n)
    {}
    virtual ~TagCAP()
    {}
    virtual void* allocate() = 0;
    virtual void deallocate(void* p) = 0;
    virtual int getType() = 0;
    bool has(MeshEntity* e)
    {
      int count = tagContainer.count(e);
      return count > 0;
    }
    void set(MeshEntity* e, void* p)
    {
      tagContainer[e] = p;
    }
    void* get(MeshEntity* e)
    {
      if ( ! has(e))
        set(e,this->allocate());
      void* p = tagContainer[e];
      return p;
    }
    void remove(MeshEntity* e)
    {
      this->deallocate(this->get(e));
      tagContainer.erase(e);
    }
    MeshDatabaseInterface* mesh;
    int count;
    std::string name;
    std::map<MeshEntity*, void*> tagContainer;
};

class DoubleTagCAP : public TagCAP
{
  public:
    DoubleTagCAP(MeshDatabaseInterface* m, const char* name, int c):
      TagCAP(m, name,c)
    {}
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
      double* internal = static_cast<double*>(TagCAP::get(e));
      for (int i=0; i < count; ++i)
        p[i] = internal[i];
    }
    void set(MeshEntity* e, double const* p)
    {
      double* internal = static_cast<double*>(TagCAP::get(e));
      for (int i=0; i < count; ++i)
        internal[i] = p[i];
    }
};

class IntTagCAP : public TagCAP
{
  public:
    IntTagCAP(MeshDatabaseInterface* m, const char* name, int c):
      TagCAP(m, name,c)
    {}
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
      int* internal = static_cast<int*>(TagCAP::get(e));
      for (int i=0; i < count; ++i)
        p[i] = internal[i];
    }
    void set(MeshEntity* e, int const* p)
    {
      int* internal = static_cast<int*>(TagCAP::get(e));
      for (int i=0; i < count; ++i)
        internal[i] = p[i];
    }
};


MeshTag* MeshCAP::createDoubleTag(const char* name, int size)
{
  TagCAP* tag = new DoubleTagCAP(meshInterface, name,size);
  tags.push_back(tag);
  return reinterpret_cast<MeshTag*>(tag);
}

MeshTag* MeshCAP::createIntTag(const char* name, int size)
{
  TagCAP* tag = new IntTagCAP(meshInterface, name,size);
  tags.push_back(tag);
  return reinterpret_cast<MeshTag*>(tag);
}

MeshTag* MeshCAP::createLongTag(const char* name, int size)
{
  (void)name;
  (void)size;
  apf::fail("MeshCAP::createLongTag called!\n");
  return 0;
}

MeshTag* MeshCAP::findTag(const char* name)
{
  for (size_t i=0; i < tags.size(); ++i)
    if (tags[i]->name == name)
      return reinterpret_cast<MeshTag*>(tags[i]);
  return 0;
}

void MeshCAP::destroyTag(MeshTag* tag)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  tags.erase(std::find(tags.begin(),tags.end(),tagCap));
  delete tagCap;;
}

void MeshCAP::renameTag(MeshTag* tag, const char*)
{
  (void)tag;
  apf::fail("MeshCAP::renameTag called!\n");
}

unsigned MeshCAP::getTagChecksum(MeshTag*,int)
{
  apf::fail("MeshCAP::getTagChecksum called!\n");
}

void MeshCAP::getTags(DynamicArray<MeshTag*>& ts)
{
  ts.setSize(tags.size());
  for (size_t i=0; i < tags.size(); ++i)
    ts[i] = reinterpret_cast<MeshTag*>(tags[i]);
}

void MeshCAP::getDoubleTag(MeshEntity* e, MeshTag* tag, double* data)
{
  DoubleTagCAP* tagCap = reinterpret_cast<DoubleTagCAP*>(tag);
  tagCap->get(e,data);
}

void MeshCAP::setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data)
{
  DoubleTagCAP* tagCap = reinterpret_cast<DoubleTagCAP*>(tag);
  tagCap->set(e,data);
}

void MeshCAP::getIntTag(MeshEntity* e, MeshTag* tag, int* data)
{
  IntTagCAP* tagCap = reinterpret_cast<IntTagCAP*>(tag);
  tagCap->get(e,data);
}

void MeshCAP::setIntTag(MeshEntity* e, MeshTag* tag, int const* data)
{
  IntTagCAP* tagCap = reinterpret_cast<IntTagCAP*>(tag);
  tagCap->set(e,data);
}

void MeshCAP::getLongTag(MeshEntity* e, MeshTag* tag, long* data)
{
  (void)e;
  (void)tag;
  (void)data;
  apf::fail("MeshCAP::getLongTag called!\n");
}

void MeshCAP::setLongTag(MeshEntity* e, MeshTag* tag, long const* data)
{
  (void)e;
  (void)tag;
  (void)data;
  apf::fail("MeshCAP::setLongTag called!\n");
}

void MeshCAP::removeTag(MeshEntity* e, MeshTag* tag)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  tagCap->remove(e);
}

bool MeshCAP::hasTag(MeshEntity* e, MeshTag* tag)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  return tagCap->has(e);
}

int MeshCAP::getTagType(MeshTag* tag)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  return tagCap->getType();
}

int MeshCAP::getTagSize(MeshTag* tag)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  return tagCap->count;
}

const char* MeshCAP::getTagName(MeshTag* tag)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  return tagCap->name.c_str();
}

bool MeshCAP::isShared(MeshEntity* e)
{
  (void)e;
  if (PCU_Comm_Peers() != 1)
    apf::fail("MeshCAP::isShared called in a parallel run!\n");
  return false;
}

bool MeshCAP::isOwned(MeshEntity* e)
{
  (void)e;
  if (PCU_Comm_Peers() != 1)
    apf::fail("MeshCAP::isOwned called in a parallel run!\n");
  return true;
}

int MeshCAP::getOwner(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::getOwner called!\n");
  return 0;
}

void MeshCAP::getRemotes(MeshEntity* e, Copies& remotes)
{
  (void)e;
  (void)remotes;
  apf::fail("MeshCAP::getRemotes called!\n");
}

void MeshCAP::getResidence(MeshEntity* e, Parts& residence)
{
  (void)e;
  (void)residence;
  apf::fail("MeshCAP::getResidence called!\n");
}

int MeshCAP::getId()
{
  apf::fail("MeshCAP::getId called!\n");
}

void MeshCAP::migrate(Migration* plan)
{
  (void)plan;
  apf::fail("MeshCAP::migrate called!\n");
}

void MeshCAP::getMatches(MeshEntity* e, Matches& m)
{
  (void)e;
  (void)m;
  apf::fail("MeshCAP::getMatches called!\n");
}

void MeshCAP::getDgCopies(MeshEntity* e, DgCopies& dgCopies, ModelEntity* me)
{
  (void)e;
  (void)dgCopies;
  (void)me;
  apf::fail("MeshCAP::getDgCopies called!\n");
}

Mesh2* createMesh(capMesh* mesh)
{
  (void)mesh;
  apf::fail("MeshCAP::createMesh called!\n");
  return 0;
}

MeshEntity* castEntity(capEntity* entity)
{
  (void)entity;
  apf::fail("MeshCAP::castEntity called!\n");
  return 0;
}

Mesh2* createMesh(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb)
{
  MeshCAP* m = new MeshCAP(mdb, gdb);
  m->init(getLagrange(1));
  return m;
}


}//namespace apf
