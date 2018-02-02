#include "apfCAP.h"
#include <apf.h>
#include <apfShape.h>
#include <gmi.h>
#include <cstdlib>
#include <pcu_util.h>
#include <algorithm>


namespace apf {


MeshCAP::MeshCAP(MeshDatabaseInterface* mdb):
  meshInterface(mdb)
{
  PCU_ALWAYS_ASSERT(meshInterface);

  std::size_t numRegions = 0;
  meshInterface->get_num_topos(TOPO_REGION, numRegions);
  d = numRegions ? 3 : 2;
  iterDim = -1;
  model = 0;
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

class CapstoneEntity
{
  public:
    CapstoneEntity(M_MTopo inEnt):
      ent(inEnt) {}
    M_MTopo get()
    {
      return ent;
    }
  private:
    M_MTopo ent;
};

MeshEntity* MeshCAP::iterate(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);
  if (meshInterface->iterator_end(*miter))
    return 0;
  meshInterface->iterator_next(*miter);
  CapstoneEntity* ce = new CapstoneEntity(meshInterface->iterator_value(*miter));
  return reinterpret_cast<MeshEntity*>(ce);
}

void MeshCAP::end(MeshIterator* it)
{
  MeshSmartIterator* miter = reinterpret_cast<MeshSmartIterator*>(it);
  delete miter;
}

bool MeshCAP::isShared(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::isShared called!\n");
  return false;
}

bool MeshCAP::isOwned(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::isOwned called!\n");
  return false;
}

int MeshCAP::getOwner(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::getOwner called!\n");
  return 0;
}

bool MeshCAP::hasAdjacency(int from_dim, int to_dim)
{
  (void)from_dim;
  (void)to_dim;
  apf::fail("MeshCAP::hasAdjacency called!\n");
  return false;
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

void MeshCAP::getAdjacent(MeshEntity* e,
    int dimension,
    DynamicArray<MeshEntity*>& adjacent)
{
  (void)e;
  (void)dimension;
  (void)adjacent;
  apf::fail("MeshCAP::getAdjacent called!\n");
}

int MeshCAP::getDownward(MeshEntity* e,
    int dimension,
    MeshEntity** down)
{
  (void)e;
  (void)dimension;
  (void)down;
  apf::fail("MeshCAP::getDownward called!\n");
  return 0;
}

int MeshCAP::countUpward(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::countUpward called!\n");
  return 0;
}

MeshEntity* MeshCAP::getUpward(MeshEntity* e, int i)
{
  (void)e;
  (void)i;
  apf::fail("MeshCAP::getUpward called!\n");
  return 0;
}

void MeshCAP::getUp(MeshEntity* e, Up& up)
{
  (void)e;
  (void)up;
  apf::fail("MeshCAP::getUp called!\n");
}

bool MeshCAP::hasUp(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::hasUp called!\n");
  return false;
}

void MeshCAP::getPoint_(MeshEntity* e, int node, Vector3& point)
{
  (void)node;
  CapstoneEntity* ce = reinterpret_cast<CapstoneEntity*>(e);
  M_MTopo topo = ce->get();
  if (meshInterface->is_vertex(topo))
    meshInterface->get_vertex_coord(topo, &(point[0]));
  else
    apf::fail("MeshCAP::getPoint_ is called for entity other than vertex!\n");
}

void MeshCAP::getParam(MeshEntity* e, Vector3& point)
{
  (void)e;
  (void)point;
  apf::fail("MeshCAP::getParam called!\n");
}

Mesh::Type MeshCAP::getType(MeshEntity* e)
{
  CapstoneEntity* ce = reinterpret_cast<CapstoneEntity*>(e);
  M_MTopo topo = ce->get();
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

MeshTag* MeshCAP::createDoubleTag(const char* name, int size)
{
  (void)name;
  (void)size;
  apf::fail("MeshCAP::createDoubleTag called!\n");
  return 0;
}

MeshTag* MeshCAP::createIntTag(const char* name, int size)
{
  (void)name;
  (void)size;
  apf::fail("MeshCAP::createIntTag called!\n");
  return 0;
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
  (void)name;
  apf::fail("MeshCAP::findTag called!\n");
  return 0;
}

void MeshCAP::destroyTag(MeshTag* tag)
{
  (void)tag;
  apf::fail("MeshCAP::destroyTag called!\n");
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
  (void)ts;
  apf::fail("MeshCAP::getTags called!\n");
}

void MeshCAP::getDoubleTag(MeshEntity* e, MeshTag* tag, double* data)
{
  (void)e;
  (void)tag;
  (void)data;
  apf::fail("MeshCAP::getDoubleTag called!\n");
}

void MeshCAP::setDoubleTag(MeshEntity* e, MeshTag* tag, double const* data)
{
  (void)e;
  (void)tag;
  (void)data;
  apf::fail("MeshCAP::setDoubleTag called!\n");
}

void MeshCAP::getIntTag(MeshEntity* e, MeshTag* tag, int* data)
{
  (void)e;
  (void)tag;
  (void)data;
  apf::fail("MeshCAP::getIntTag called!\n");
}

void MeshCAP::setIntTag(MeshEntity* e, MeshTag* tag, int const* data)
{
  (void)e;
  (void)tag;
  (void)data;
  apf::fail("MeshCAP::setIntTag called!\n");
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
  (void)e;
  (void)tag;
  apf::fail("MeshCAP::removeTag called!\n");
}

bool MeshCAP::hasTag(MeshEntity* e, MeshTag* tag)
{
  (void)e;
  (void)tag;
  apf::fail("MeshCAP::hasTag called!\n");
}

int MeshCAP::getTagType(MeshTag* tag)
{
  (void)tag;
  apf::fail("MeshCAP::getTagType called!\n");
}

int MeshCAP::getTagSize(MeshTag* tag)
{
  (void)tag;
  apf::fail("MeshCAP::getTagSize called!\n");
}

const char* MeshCAP::getTagName(MeshTag* tag)
{
  (void)tag;
  apf::fail("MeshCAP::getTagName called!\n");
}

ModelEntity* MeshCAP::toModel(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::toModel called!\n");
}

gmi_model* MeshCAP::getModel()
{
  apf::fail("MeshCAP::getModel called!\n");
}

void MeshCAP::migrate(Migration* plan)
{
  (void)plan;
  apf::fail("MeshCAP::migrate called!\n");
}

int MeshCAP::getId()
{
  apf::fail("MeshCAP::getId called!\n");
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

void MeshCAP::verify()
{
  apf::fail("MeshCAP::verify called!\n");
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

void MeshCAP::setPoint_(MeshEntity * me, int node, Vector3 const & p)
{
  (void)me;
  (void)node;
  (void)p;
  apf::fail("MeshCAP::setPoint_ called!\n");
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

Mesh2* createMesh(MeshDatabaseInterface* mdb)
{
  MeshCAP* m = new MeshCAP(mdb);
  m->init(getLagrange(1));
  return m;
}


}//namespace apf
