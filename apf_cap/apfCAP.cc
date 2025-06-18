#include "apfCAP.h"

#include <algorithm>
#include <cstdlib>

#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <mth.h>
#include <pcu_util.h>
#include <PCU.h>

#include <CreateMG_Framework_Mesh.h>
#ifdef HAVE_CAPSTONE_SIZINGMETRICTOOL
#include <CreateMG_SizingMetricTool.h>
#endif

using namespace CreateMG;
namespace MeshMG = CreateMG::Mesh;

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

bool hasCAP() noexcept {
  return true;
}

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

class TagCAP;

class MeshCAP : public Mesh2
{
  public:
    MeshCAP(MDBI* mdb, GDBI* gdb);
    MeshCAP(gmi_model* mdl, MDBIP mdb);
    virtual ~MeshCAP();
    /* --------------------------------------------------------------------- */
    /* Category 00: General Mesh APIs */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    int getDimension();
    std::size_t count(int dimension);
    Type getType(MeshEntity* e);
    void verify();
    // OPTIONAL Member Functions //
    void writeNative(const char* fileName);
    void destroyNative();

    /* --------------------------------------------------------------------- */
    /* Category 01: General Getters and Setters for vertex coordinates */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    void getPoint_(MeshEntity* e, int, Vector3& point);
    void setPoint_(MeshEntity* e, int, Vector3 const& p);
    void getParam(MeshEntity* e, Vector3& p);
    void setParam(MeshEntity* e, Vector3 const& point);

    /* --------------------------------------------------------------------- */
    /* Category 02: Iterators */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    MeshIterator* begin(int dimension);
    MeshEntity* iterate(MeshIterator* it);
    void end(MeshIterator* it);
    void increment(MeshIterator* it);
    bool isDone(MeshIterator* it);
    MeshEntity* deref(MeshIterator* it);

    /* --------------------------------------------------------------------- */
    /* Category 03: Adjacencies */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    void getAdjacent(MeshEntity* e, int dimension, Adjacent& adjacent);
    int getDownward(MeshEntity* e, int dimension, MeshEntity** adjacent);
    MeshEntity* getUpward(MeshEntity* e, int i);
    bool hasUp(MeshEntity* e);
    // OPTIONAL Member Functions //
    bool hasAdjacency(int from_dim, int to_dim);
    void createAdjacency(int from_dim, int to_dim);
    void deleteAdjacency(int from_dim, int to_dim);
    void getUp(MeshEntity* e, Up& up);
    int countUpward(MeshEntity* e);

    /* --------------------------------------------------------------------- */
    /* Category 04: CAD model inquires */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    ModelEntity* toModel(MeshEntity* e);
    // OPTIONAL Member Functions //
    gmi_model* getModel();
    void setModel(gmi_model* newModel);
    void setModelEntity(MeshEntity* e, ModelEntity* me);

    /* --------------------------------------------------------------------- */
    /* Category 05: Entity Creation/Deletion */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    MeshEntity* createVert_(ModelEntity* me);
    MeshEntity* createEntity_(int type, ModelEntity* me, MeshEntity** down);
    void destroy_(MeshEntity* e);

    /* --------------------------------------------------------------------- */
    /* Category 06: Attachable Data Functionality */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    MeshTag* createDoubleTag(const char* name, int size);
    MeshTag* createIntTag(const char* name, int size);
    MeshTag* createLongTag(const char* name, int size);
    MeshTag* findTag(const char* name);
    void destroyTag(MeshTag* t);
    void renameTag(MeshTag* t, const char* name);
    void getTags(DynamicArray<MeshTag*>& tags);
    /* void getTag(MeshEntity* e, MeshTag* t, void* data); */
    /* void setTag(MeshEntity* e, MeshTag* t, void const* data); */
    void getDoubleTag(MeshEntity* e, MeshTag* tag, double* data);
    void setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data);
    void getIntTag(MeshEntity* e, MeshTag* tag, int* data);
    void setIntTag(MeshEntity* e, MeshTag* tag, int const* data);
    void getLongTag(MeshEntity* e, MeshTag* tag, long* data);
    void setLongTag(MeshEntity* e, MeshTag* tag, long const* data);
    void removeTag(MeshEntity* e, MeshTag* t);
    bool hasTag(MeshEntity* e, MeshTag* t);
    int getTagType(MeshTag* t);
    int getTagSize(MeshTag* t);
    const char* getTagName(MeshTag* t);
    unsigned getTagChecksum(MeshTag*,int);


    /* --------------------------------------------------------------------- */
    /* Category 07: Distributed Meshes */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    bool isShared(MeshEntity* e);
    bool isGhost(MeshEntity*) { return false; }
    bool isGhosted(MeshEntity*) { return false; }
    bool isOwned(MeshEntity* e);
    int getOwner(MeshEntity* e);
    void getRemotes(MeshEntity* e, Copies& remotes);
    void getResidence(MeshEntity* e, Parts& residence);
    int getId();
    void setResidence(MeshEntity*, Parts&) {}
    void acceptChanges() {}
    // OPTIONAL Member Functions //
    void deleteGhost(MeshEntity*) {}
    void addGhost(MeshEntity*, int, MeshEntity*) {}
    int getGhosts(MeshEntity*, Copies&) { return 0; }
    void migrate(Migration* plan);
    void setRemotes(MeshEntity*, Copies&) {}
    void addRemote(MeshEntity*, int, MeshEntity*) {}
    void clearRemotes(MeshEntity*) {}


    /* --------------------------------------------------------------------- */
    /* Category 08: Periodic Meshes */
    /* --------------------------------------------------------------------- */
    // REQUIRED Member Functions //
    bool hasMatching() { return false; }
    void getMatches(MeshEntity* e, Matches& m);
    // OPTIONAL Member Functions //
    void addMatch(MeshEntity*, int, MeshEntity*) {}
    void clearMatches(MeshEntity*) {}
    void clear_() {}
    void getDgCopies(MeshEntity* e, DgCopies& dgCopies, ModelEntity* me);

    MDBI* getMesh() { return meshInterface; }
  protected:
    gmi_model* model;
    MDBIP meshOwner;
    MDBI* meshInterface;
    int d;
    std::vector<TagCAP*> tags;
};

MeshCAP::MeshCAP(MDBI* mdb, GDBI* gdb): meshInterface(mdb) {
  model = gmi_import_cap(gdb);
  PCU_ALWAYS_ASSERT(mdb);
  MG_API_CALL(mdb, get_dimension(d));
}

MeshCAP::MeshCAP(gmi_model* mdl, MDBIP mdb):
  model(mdl), meshOwner(mdb), meshInterface(mdb.get())
{
  PCU_ALWAYS_ASSERT(mdl);
  PCU_ALWAYS_ASSERT(mdb.is_valid());
  MG_API_CALL(mdb.get(), get_dimension(d));
}

MeshCAP::~MeshCAP()
{
  if (model) gmi_destroy(model);
  model = nullptr;
  meshInterface = nullptr;
}

int MeshCAP::getDimension()
{
  return d;
}

std::size_t MeshCAP::count(int dimension)
{
  std::size_t count = 0;
  if (dimension == 0)
    meshInterface->get_num_topos(MeshMG::TOPO_VERTEX, count);
  if (dimension == 1)
    meshInterface->get_num_topos(MeshMG::TOPO_EDGE, count);
  if (dimension == 2)
    meshInterface->get_num_topos(MeshMG::TOPO_FACE, count);
  if (dimension == 3)
    meshInterface->get_num_topos(MeshMG::TOPO_REGION, count);
  return count;
}

Mesh::Type MeshCAP::getType(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  MeshMG::MeshShape topoShape;
  meshInterface->get_topo_shape(topo, topoShape);
  if (topoShape == MeshMG::SHAPE_NODE) return VERTEX;
  else if (topoShape == MeshMG::SHAPE_SEGMENT) return EDGE;
  else if (topoShape == MeshMG::SHAPE_TRIANGLE) return TRIANGLE;
  else if (topoShape == MeshMG::SHAPE_QUAD) return QUAD;
  else if (topoShape == MeshMG::SHAPE_TETRA) return TET;
  else if (topoShape == MeshMG::SHAPE_HEX) return HEX;
  else if (topoShape == MeshMG::SHAPE_PRISM) return PRISM;
  else if (topoShape == MeshMG::SHAPE_PYRAMID) return PYRAMID;
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
  MeshMG::GeometryTopoType gtype;
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
  auto miter = new MeshMG::MeshSmartIterator(meshInterface);
  if (dimension == 0)
    meshInterface->get_topo_iterator(MeshMG::TOPO_VERTEX, *miter);
  if (dimension == 1)
    meshInterface->get_topo_iterator(MeshMG::TOPO_EDGE, *miter);
  if (dimension == 2)
    meshInterface->get_topo_iterator(MeshMG::TOPO_FACE, *miter);
  if (dimension == 3)
    meshInterface->get_topo_iterator(MeshMG::TOPO_REGION, *miter);
  meshInterface->iterator_begin(*miter);
  return reinterpret_cast<MeshIterator*>(miter);
}

/* NOTE: miter is located at the first item in the list, therefore
 * iterate has to return it before calling iterator_next on miter
 */
MeshEntity* MeshCAP::iterate(MeshIterator* it)
{
  auto miter = reinterpret_cast<MeshMG::MeshSmartIterator*>(it);

  M_MTopo topo = meshInterface->iterator_value(*miter);

  if (!meshInterface->iterator_end(*miter))
    meshInterface->iterator_next(*miter);
  else
    return 0;

  return toEntity(topo);
}

void MeshCAP::end(MeshIterator* it)
{
  auto miter = reinterpret_cast<MeshMG::MeshSmartIterator*>(it);
  delete miter;
}

void MeshCAP::increment(MeshIterator* it)
{
  auto miter = reinterpret_cast<MeshMG::MeshSmartIterator*>(it);
  meshInterface->iterator_next(*miter);
  it = reinterpret_cast<MeshIterator*>(miter);
}

bool MeshCAP::isDone(MeshIterator* it)
{
  auto miter = reinterpret_cast<MeshMG::MeshSmartIterator*>(it);
  return meshInterface->iterator_end(*miter);
}

MeshEntity* MeshCAP::deref(MeshIterator* it)
{
  auto miter = reinterpret_cast<MeshMG::MeshSmartIterator*>(it);
  M_MTopo topo = meshInterface->iterator_value(*miter);
  return toEntity(topo);
}

void MeshCAP::getAdjacent(MeshEntity* e,
    int dimension,
    DynamicArray<MeshEntity*>& adjacent)
{
  M_MTopo topo = fromEntity(e);
  MeshMG::MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<M_MTopo> adjTopos;
  if (apf::getDimension(this, e) == dimension)
  {
    adjacent.setSize(1);
    adjacent[0] = e;
    return;
  }
  if (type == MeshMG::TOPO_VERTEX)
  {
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_EDGE, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_FACE, adjTopos);
    if (dimension == 3)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_REGION, adjTopos);
  }
  if (type == MeshMG::TOPO_EDGE)
  {
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_VERTEX, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_FACE, adjTopos);
    if (dimension == 3)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_REGION, adjTopos);
  }
  if (type == MeshMG::TOPO_FACE)
  {
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_VERTEX, adjTopos);
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_EDGE, adjTopos);
    if (dimension == 3)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_REGION, adjTopos);
  }
  if (type == MeshMG::TOPO_REGION)
  {
    if (dimension == 0)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_VERTEX, adjTopos);
    if (dimension == 1)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_EDGE, adjTopos);
    if (dimension == 2)
      meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_FACE, adjTopos);
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
  MeshMG::MeshTopo type;
  meshInterface->get_topo_type(topo, type);

  int from;
  if (type == MeshMG::TOPO_VERTEX)
    from = 0;
  else if(type == MeshMG::TOPO_EDGE)
    from = 1;
  else if(type == MeshMG::TOPO_FACE)
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
    meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_VERTEX, adjTopos);
  if (dimension == 1)
    meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_EDGE, adjTopos);
  if (dimension == 2)
    meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_FACE, adjTopos);

  for (std::size_t i = 0; i < adjTopos.size(); i++)
    down[i] = toEntity(adjTopos[i]);

  PCU_ALWAYS_ASSERT(adjTopos.size() == (std::size_t)Mesh::adjacentCount[getType(e)][dimension]);
  return adjTopos.size();
}

MeshEntity* MeshCAP::getUpward(MeshEntity* e, int i)
{
  M_MTopo topo = fromEntity(e);
  MeshMG::MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<std::size_t> adjId;
  if (type == MeshMG::TOPO_VERTEX)
  {
    meshInterface->get_adjacency_id_vector(topo, MeshMG::TOPO_EDGE, adjId);
    M_MTopo topo = meshInterface->get_topo_by_id(MeshMG::TOPO_EDGE, adjId[i]);
    return toEntity(topo);
  }
  if (type == MeshMG::TOPO_EDGE)
  {
    meshInterface->get_adjacency_id_vector(topo, MeshMG::TOPO_FACE, adjId);
    M_MTopo topo = meshInterface->get_topo_by_id(MeshMG::TOPO_FACE, adjId[i]);
    return toEntity(topo);
  }
  if (type == MeshMG::TOPO_FACE)
  {
    meshInterface->get_adjacency_id_vector(topo, MeshMG::TOPO_REGION, adjId);
    M_MTopo topo = meshInterface->get_topo_by_id(MeshMG::TOPO_REGION, adjId[i]);
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
  MeshMG::MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<M_MTopo> adjTopos;
  if (type == MeshMG::TOPO_VERTEX)
    meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_EDGE, adjTopos);
  if (type == MeshMG::TOPO_EDGE)
    meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_FACE, adjTopos);
  if (type == MeshMG::TOPO_FACE)
    meshInterface->get_adjacency_vector(topo, MeshMG::TOPO_REGION, adjTopos);
  up.n = adjTopos.size();
  for (int i = 0; i < up.n; i++) {
    up.e[i] = toEntity(adjTopos[i]);
  }
}

int MeshCAP::countUpward(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  MeshMG::MeshTopo type;
  meshInterface->get_topo_type(topo, type);
  std::vector<std::size_t> adjTopos;
  if (type == MeshMG::TOPO_VERTEX)
    meshInterface->get_adjacency_id_vector(topo, MeshMG::TOPO_EDGE, adjTopos);
  if (type == MeshMG::TOPO_EDGE)
    meshInterface->get_adjacency_id_vector(topo, MeshMG::TOPO_FACE, adjTopos);
  if (type == MeshMG::TOPO_FACE)
    meshInterface->get_adjacency_id_vector(topo, MeshMG::TOPO_REGION, adjTopos);
  return (int)adjTopos.size();
}

ModelEntity* MeshCAP::toModel(MeshEntity* e)
{
  M_MTopo topo = fromEntity(e);
  M_GTopo gtopo;
  MeshMG::GeometryTopoType gtype;
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

static MeshMG::GeometryTopoType getCapGeomType(int d)
{
  MeshMG::GeometryTopoType gtype = MeshMG::GVERTEX;
  switch (d) {
    case 0:
      gtype = MeshMG::GVERTEX;
      break;
    case 1:
      gtype = MeshMG::GEDGE;
      break;
    case 2:
      gtype = MeshMG::GFACE;
      break;
    case 3:
      gtype = MeshMG::GREGION;
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
  MeshMG::GeometryTopoType gtype = getCapGeomType(d);
  meshInterface->create_vertex(xyz, vertex, gtype, gtopo);
  return toEntity(vertex);
}

static MeshMG::MeshShape getCapShape(int type)
{
  MeshMG::MeshShape shape = MeshMG::SHAPE_UNKNOWN;
  switch (type) {
    case Mesh::VERTEX:
      shape = MeshMG::SHAPE_NODE;
      break;
    case Mesh::EDGE:
      shape = MeshMG::SHAPE_SEGMENT;
      break;
    case Mesh::TRIANGLE:
      shape = MeshMG::SHAPE_TRIANGLE;
      break;
    case Mesh::QUAD:
      shape = MeshMG::SHAPE_QUAD;
      break;
    case Mesh::TET:
      shape = MeshMG::SHAPE_TETRA;
      break;
    case Mesh::HEX:
      shape = MeshMG::SHAPE_HEX;
      break;
    case Mesh::PRISM:
      shape = MeshMG::SHAPE_PRISM;
      break;
    case Mesh::PYRAMID:
      shape = MeshMG::SHAPE_PYRAMID;
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
  MeshMG::MeshShape shape = getCapShape(type);
  if ( !me ) {
    meshInterface->create_topo(shape, mtopos, topo);
  }
  // if model is not 0 figure out its type
  else {
    int d = getModelType(me);
    MeshMG::GeometryTopoType gtype = getCapGeomType(d);
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
    TagCAP(MDBI* m, const char* n, int c):
      mesh(m),
      count(c),
      name(n)
    {}
    virtual ~TagCAP() {}
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
    MDBI* mesh;
    int count;
    std::string name;
    std::map<MeshEntity*, void*> tagContainer;
};

class DoubleTagCAP : public TagCAP
{
  public:
    DoubleTagCAP(MDBI* m, const char* name, int c):
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
    IntTagCAP(MDBI* m, const char* name, int c):
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

Mesh2* createCapMesh(MDBI* mdb, GDBI* gdb, pcu::PCU *PCUObj) {
  MeshCAP* m = new MeshCAP(mdb, gdb);
  m->init(getLagrange(1), PCUObj);
  return m;
}

Mesh2* createCapMesh(
  gmi_model* model, const char* meshname, pcu::PCU* PCUObj
) {
  GDBI* gdb = gmi_export_cap(model);
  AppContext* ctx = gdb->get_context();
  MDBIP mdp = get_context_mesh_database_interface(ctx);
  M_MModel mmodel;
  MG_API_CALL(mdp.get(), get_model_by_name(meshname, mmodel));
  MG_API_CALL(mdp.get(), set_current_model(mmodel));
  MeshCAP* m = new MeshCAP(model, mdp);
  m->init(getLagrange(1), PCUObj);
  return m;
}

bool has_smoothCAPAnisoSizes(void) noexcept {
#ifdef HAVE_CAPSTONE_SIZINGMETRICTOOL
  return true;
#else
  return false;
#endif
}

bool smoothCAPAnisoSizes(apf::Mesh2* mesh, std::string analysis,
  apf::Field* scales, apf::Field* frames) {
#ifdef HAVE_CAPSTONE_SIZINGMETRICTOOL
  // Ensure input is a MeshCAP.
  apf::MeshCAP* m = dynamic_cast<apf::MeshCAP*>(mesh);
  if (!m) {
    lion_eprint(1, "ERROR: smoothCAPAnisoSizes: mesh is not an apf::MeshCAP*\n");
    return false;
  }

  // Extract metric tensors from MeshAdapt frames and scales.
  std::vector<Metric6> sizing6(m->count(0));
  apf::Matrix3x3 Q;
  apf::Vector3 H;
  apf::MeshIterator* it = m->begin(0);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    apf::getVector(scales, e, 0, H); // Desired element lengths.
    apf::getMatrix(frames, e, 0, Q); // MeshAdapt uses column vectors.
    apf::Matrix3x3 L(1.0/(H[0]*H[0]), 0, 0,
      0, 1.0/(H[1]*H[1]), 0,
      0, 0, 1.0/(H[2]*H[2]));
    apf::Matrix3x3 t = Q * L * apf::transpose(Q); // Invert orthogonal frames.
    size_t id;
    MG_API_CALL(m->getMesh(), get_id(fromEntity(e), id));
    PCU_DEBUG_ASSERT(id != 0);
    --id;
    sizing6[id][0] = t[0][0];
    sizing6[id][1] = t[0][1];
    sizing6[id][2] = t[0][2];
    sizing6[id][3] = t[1][1];
    sizing6[id][4] = t[1][2];
    sizing6[id][5] = t[2][2];
  }
  m->end(it);
  auto smooth_tool = get_sizing_metric_tool(m->getMesh()->get_context(),
    "CreateSmoothingBase");
  if (smooth_tool == nullptr) {
    lion_eprint(1, "ERROR: Unable to find \"CreateSmoothingBase\"\n");
    return false;
  }
  smooth_tool->set_context(m->getMesh()->get_context());
  M_MModel mmodel;
  MG_API_CALL(m->getMesh(), get_current_model(mmodel));
  smooth_tool->set_metric(mmodel, "sizing6", sizing6);
  std::vector<Metric6> ometric;
  smooth_tool->smooth_metric(mmodel, analysis, "sizing6", ometric);
  it = m->begin(0);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    size_t id;
    MG_API_CALL(m->getMesh(), get_id(fromEntity(e), id));
    PCU_DEBUG_ASSERT(id != 0);
    --id;
    const Metric6& m = ometric[id];
    apf::Matrix3x3 t(m[0], m[1], m[2],
      m[1], m[3], m[4],
      m[2], m[4], m[5]);
    int n = apf::eigen(t, &Q[0], &H[0]); // Eigenvectors in rows of Q.
    PCU_DEBUG_ASSERT(n == 3);
    Q = apf::transpose(Q); // Put eigenvectors back into columns for MeshAdapt.
    for (int i = 0; i < 3; ++i) {
      H[i] = 1.0/sqrt(H[i]);
    }
    apf::setMatrix(frames, e, 0, Q);
    apf::setVector(scales, e, 0, H);
  }
  m->end(it);
  return true;
#else
  (void) mesh;
  (void) analysis;
  (void) scales;
  (void) frames;
  apf::fail("smoothCAPAnisoSizes: Capstone does not have SizingMetricTool.");
#endif
}

}//namespace apf
