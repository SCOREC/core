#include "apfCAP.h"
#include <apf.h>
#include <apfShape.h>
#include <gmi.h>
#include <cstdlib>
#include <pcu_util.h>
#include <algorithm>

namespace apf {


MeshCAP::MeshCAP(capMesh* m):
  mesh(m)
{
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
  (void)dimension;
  apf::fail("MeshCAP::count called!\n");
  return 0;
}

MeshIterator* MeshCAP::begin(int dimension)
{
  (void)dimension;
  apf::fail("MeshCAP::begin called!\n");
  return 0;
}

MeshEntity* MeshCAP::iterate(MeshIterator* it)
{
  (void)it;
  apf::fail("MeshCAP::iterate called!\n");
  return 0;
}

void MeshCAP::end(MeshIterator* it)
{
  (void)it;
  apf::fail("MeshCAP::end called!\n");
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
  (void)e;
  (void)node;
  (void)point;
  apf::fail("MeshCAP::getPoint_ called!\n");
}

void MeshCAP::getParam(MeshEntity* e, Vector3& point)
{
  (void)e;
  (void)point;
  apf::fail("MeshCAP::getParam called!\n");
}

Mesh::Type MeshCAP::getType(MeshEntity* e)
{
  (void)e;
  apf::fail("MeshCAP::getType called!\n");
  return Mesh::VERTEX;
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

}//namespace apf
