#include "apfCAP.h"
#include <apf.h>
#include <apfShape.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <cstdlib>
#include <pcu_util.h>
#include <algorithm>


namespace apf {

/* These tables are taken from ../mds/mds.c
 * Here they are used to figure out the order
 * of verts in a tet given the downward faces,
 * for example.
 */
static int const t01[] = {0,1,1,2,2,0};
static int const t10[] = {2,0,0,1,1,2};

static int const q01[] = {};
static int const q10[] = {};

static int const T01[] = {0,1,1,2,2,0
                         ,0,3,1,3,2,3};
static int const T10[] = {2,0,0,1,1,2,3,4};
static int const T12[] = {0,1,2
                         ,0,4,3
                         ,1,5,4
                         ,2,3,5};
static int const T21[] = {0,1
                         ,0,2
                         ,0,3
                         ,1,3
                         ,1,2
                         ,2,3};

static int const* convs[Mesh::TYPES][4][4] =
{{{0,0  ,0,0},{0  ,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,0  ,0,0},{0  ,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,t01,0,0},{t10,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,q01,0,0},{q10,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,T01,0,0},{T10,0,T12,0},{0,T21,0,0},{0,0,0,0}}
,{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}
,{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}
,{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}
};

// entity type-dimension count table
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

MeshCAP::MeshCAP(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb):
  meshInterface(mdb), geomInterface(gdb)
{
  PCU_ALWAYS_ASSERT(meshInterface);

  std::size_t numRegions = 0;
  meshInterface->get_num_topos(TOPO_REGION, numRegions);
  d = numRegions ? 3 : 2;
  iterDim = -1;
  model = gmi_import_cap(geomInterface);
}

MeshCAP::~MeshCAP()
{
  if (model) destroyNative();
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
  apf::verify(this);
  // other verification using Capstone APIs can go here.
}

void MeshCAP::writeNative(const char* fileName)
{
  (void)fileName;
  apf::fail("MeshCAP::writeNative called!\n");
}

void MeshCAP::destroyNative()
{
  if (model) gmi_destroy(model);
  model = nullptr;
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
  double u, v;
  GeometryTopoType gtype;
  meshInterface->get_vertex_uv_parameters(topo, u, v, gtype);
  point = Vector3(u, v, 0.);
}

void MeshCAP::setParam(MeshEntity* e, Vector3 const& point)
{
  M_MTopo topo = fromEntity(e);
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

  int from;
  if (type == TOPO_VERTEX)
    from = 0;
  else if(type == TOPO_EDGE)
    from = 1;
  else if(type == TOPO_FACE)
    from = 2;
  else
    from = 3;

  std::vector<M_MTopo> adjTopos;
  if (from == dimension)
  {
    down[0] = e;
    return 1;
  }

  PCU_ALWAYS_ASSERT(from > dimension);
  if (dimension == 0)
    meshInterface->get_adjacency_vector(topo, TOPO_VERTEX, adjTopos);
  if (dimension == 1)
    meshInterface->get_adjacency_vector(topo, TOPO_EDGE, adjTopos);
  if (dimension == 2)
    meshInterface->get_adjacency_vector(topo, TOPO_FACE, adjTopos);

  for (std::size_t i = 0; i < adjTopos.size(); i++)
    down[i] = toEntity(adjTopos[i]);

  PCU_ALWAYS_ASSERT(adjTopos.size() == (std::size_t)Mesh::adjacentCount[getType(e)][dimension]);
  return adjTopos.size();
}

MeshEntity* MeshCAP::getUpward(MeshEntity* e, int i)
{
  M_MTopo topo = fromEntity(e);
  MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<std::size_t> adjId;
  if (type == TOPO_VERTEX)
  {
    meshInterface->get_adjacency_id_vector(topo, TOPO_EDGE, adjId);
    M_MTopo topo = meshInterface->get_topo_by_id(TOPO_EDGE, adjId[i]);
    return toEntity(topo);
  }
  if (type == TOPO_EDGE)
  {
    meshInterface->get_adjacency_id_vector(topo, TOPO_FACE, adjId);
    M_MTopo topo = meshInterface->get_topo_by_id(TOPO_FACE, adjId[i]);
    return toEntity(topo);
  }
  if (type == TOPO_FACE)
  {
    meshInterface->get_adjacency_id_vector(topo, TOPO_REGION, adjId);
    M_MTopo topo = meshInterface->get_topo_by_id(TOPO_REGION, adjId[i]);
    return toEntity(topo);
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
  up.n = adjTopos.size();
  for (int i = 0; i < up.n; i++) {
    up.e[i] = toEntity(adjTopos[i]);
  }
}

int MeshCAP::countUpward(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<std::size_t> adjTopos;
  if (type == TOPO_VERTEX)
    meshInterface->get_adjacency_id_vector(topo, TOPO_EDGE, adjTopos);
  if (type == TOPO_EDGE)
    meshInterface->get_adjacency_id_vector(topo, TOPO_FACE, adjTopos);
  if (type == TOPO_FACE)
    meshInterface->get_adjacency_id_vector(topo, TOPO_REGION, adjTopos);
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

void MeshCAP::setModel(gmi_model* newModel)
{
  model = newModel;
}

void MeshCAP::setModelEntity(MeshEntity* e, ModelEntity* me)
{
  (void)e;
  (void)me;
  apf::fail("MeshCAP::setModelEntity called!\n");
}

static GeometryTopoType getCapGeomType(int d)
{
  GeometryTopoType gtype = GVERTEX;
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
  return gtype;
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
  int d = getModelType(me);
  GeometryTopoType gtype = getCapGeomType(d);
  meshInterface->create_vertex(xyz, vertex, gtype, gtopo);
  return toEntity(vertex);
}

static MeshShape getCapShape(int type)
{
  MeshShape shape = SHAPE_UNKNOWN;
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
  return shape;
}

static MeshEntity* commonDown(Mesh2* m, MeshEntity* a, MeshEntity* b, int dim)
{
  MeshEntity* aDown[12];
  MeshEntity* bDown[12];
  int na = m->getDownward(a, dim, aDown);
  int nb = m->getDownward(b, dim, bDown);
  PCU_ALWAYS_ASSERT( na == degree[m->getType(a)][dim] );
  PCU_ALWAYS_ASSERT( na == nb );
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++)
      if (aDown[i] == bDown[j])
      	return aDown[i];
  return 0;
}

static void stepDown(Mesh2* m, int type, int fromDim, MeshEntity** from, MeshEntity** to)
{
  PCU_ALWAYS_ASSERT(fromDim > 0);
  int toDim = fromDim - 1;
  int const* conversion = convs[type][fromDim][toDim];
  int toCount = degree[type][toDim];
  for (int i = 0; i < toCount; i++) {
    MeshEntity* a = from[conversion[2*i]];
    MeshEntity* b = from[conversion[2*i+1]];
    to[i] = commonDown(m, a, b, toDim);
  }
}

MeshEntity* MeshCAP::createEntity_(int type, ModelEntity* me, MeshEntity** down)
{
  int downType = getType(down[0]);
  std::vector<M_MTopo> mtopos;
  if (type == Mesh::TET && downType == Mesh::TRIANGLE) {
    // going from faces to verts
    // 1- go from faces to edges
    MeshEntity* downEdges[6];
    stepDown(this, Mesh::TET, 2, down, downEdges);
    // 2- go from edges to verts
    MeshEntity* downVerts[4];
    stepDown(this, Mesh::TET, 1, downEdges, downVerts);
    for (int i = 0; i < 4; i++)
      mtopos.push_back(fromEntity(downVerts[i]));
  }
  else if (type == Mesh::TRIANGLE && downType == Mesh::EDGE) {
    // going from edges to verts
    MeshEntity* downVerts[3];
    stepDown(this, Mesh::TRIANGLE, 1, down, downVerts);
    for (int i = 0; i < 3; i++)
      mtopos.push_back(fromEntity(downVerts[i]));
  }
  else if (type == Mesh::EDGE && downType == Mesh::VERTEX) {
    for (int i = 0; i < 2; i++)
      mtopos.push_back(fromEntity(down[i]));
  }
  else {
    apf::fail("MeshCAP::createEntity_ called for unsupported types!\n");
  }

  M_MTopo topo;
  MeshShape shape = getCapShape(type);
  if ( !me ) {
    meshInterface->create_topo(shape, mtopos, topo);
  }
  // if model is not 0 figure out its type
  else {
    int d = getModelType(me);
    GeometryTopoType gtype = getCapGeomType(d);
    gmi_ent* g = reinterpret_cast<gmi_ent*>(me);
    M_GTopo gtopo = fromGmiEntity(g);
    meshInterface->create_topo(shape, mtopos, topo, gtype, gtopo);
  }

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
    void rename(const char* n)
    {
      this->name = n;
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

void MeshCAP::renameTag(MeshTag* tag, const char* name)
{
  TagCAP* tagCap = reinterpret_cast<TagCAP*>(tag);
  tagCap->rename(name);
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
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::isShared called in a parallel run!\n");
  return false;
}

bool MeshCAP::isOwned(MeshEntity* e)
{
  (void)e;
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::isOwned called in a parallel run!\n");
  return true;
}

int MeshCAP::getOwner(MeshEntity* e)
{
  (void)e;
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::getOwner called in a parallel run!\n");
  return 0;
}

void MeshCAP::getRemotes(MeshEntity* e, Copies& remotes)
{
  (void)e;
  (void)remotes;
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::getRemotes called in a parallel run!\n");
}

void MeshCAP::getResidence(MeshEntity* e, Parts& residence)
{
  (void)e;
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::getResidence called in a parallel run!\n");
  residence.insert(0);
}

int MeshCAP::getId()
{
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::getId called in a parallel run!\n");
  return getPCU()->Self();
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
  if (getPCU()->Peers() != 1)
    apf::fail("MeshCAP::getMatches called in a parallel run!\n");
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

Mesh2* createMesh(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb, pcu::PCU *PCUObj)
{
  MeshCAP* m = new MeshCAP(mdb, gdb);
  m->init(getLagrange(1), PCUObj);
  return m;
}


}//namespace apf
