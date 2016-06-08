#include "apfOmega_h.h"

#include <cassert>

#include <apfMesh2.h>
#include <apfNumbering.h>
#include <PCU.h>
#include <apf.h>

namespace osh {

static void connectivityFromAPF(
    apf::Mesh* am,
    osh_t om,
    apf::Numbering* local,
    unsigned ent_dim)
{
  unsigned verts_per_elem = ent_dim + 1;
  unsigned nents = (unsigned) am->count(ent_dim);
  unsigned* conn = osh_build_ents(om, ent_dim, nents);
  apf::MeshEntity* e;
  apf::MeshIterator* it = am->begin(ent_dim);
  unsigned i = 0;
  while ((e = am->iterate(it))) {
    apf::NewArray<int> e2v;
    apf::getElementNumbers(local, e, e2v);
    for (unsigned j = 0; j < verts_per_elem; ++j)
      conn[i * verts_per_elem + j] = (unsigned) e2v[j];
    ++i;
  }
  am->end(it);
}

static void classificationFromAPF(
    apf::Mesh* am,
    osh_t om,
    unsigned ent_dim)
{
  unsigned* class_dim = osh_new_label(om, ent_dim, "class_dim", 1);
  unsigned* class_id = osh_new_label(om, ent_dim, "class_id", 1);
  apf::MeshEntity* e;
  apf::MeshIterator* it = am->begin(ent_dim);
  unsigned i = 0;
  while ((e = am->iterate(it))) {
    apf::ModelEntity* ge = am->toModel(e);
    class_dim[i] = am->getModelType(ge);
    class_id[i] = am->getModelTag(ge);
    ++i;
  }
  am->end(it);
}

static void globalFromAPF(
    apf::Mesh* am,
    osh_t om,
    unsigned ent_dim)
{
  apf::GlobalNumbering* glob_n = apf::makeGlobal(
      apf::numberOwnedDimension(am, "osh_global", ent_dim));
  apf::synchronize(glob_n);
  unsigned long* global = osh_new_global(om, ent_dim);
  apf::MeshEntity* e;
  apf::MeshIterator* it = am->begin(ent_dim);
  unsigned i = 0;
  while ((e = am->iterate(it))) {
    global[i] = apf::getNumber(glob_n, apf::Node(e, 0));
    ++i;
  }
  am->end(it);
  apf::destroyGlobalNumbering(glob_n);
}

static void fieldFromAPF(
    apf::Mesh* am,
    osh_t om,
    apf::Field* f)
{
  char const* name = apf::getName(f);
  int nc = apf::countComponents(f);
  osh_new_field(om, 0, name, nc, OSH_TRANSFER_NOT);
  double* data = osh_get_field(om, 0, name);
  apf::MeshEntity* v;
  apf::MeshIterator* it = am->begin(0);
  unsigned i = 0;
  while ((v = am->iterate(it))) {
    apf::getComponents(f, v, 0, data + i * nc);
    ++i;
  }
  am->end(it);
}

static void fieldToAPF(
    apf::Mesh* am,
    osh_t om,
    char const* name)
{
  unsigned nc = osh_components(om, 0, name);
  double* data = osh_get_field(om, 0, name);
  apf::Field* f = apf::createPackedField(am, name, nc);
  apf::MeshEntity* v;
  apf::MeshIterator* it = am->begin(0);
  unsigned i = 0;
  while ((v = am->iterate(it))) {
    apf::setComponents(f, v, 0, data + i * nc);
    ++i;
  }
  am->end(it);
}

static void fieldsFromAPF(
    apf::Mesh* am,
    osh_t om)
{
  for (int i = 0; i < am->countFields(); ++i)
    if (apf::getShape(am->getField(i)) == am->getShape())
      fieldFromAPF(am, om, am->getField(i));
}

static void fieldsToAPF(
    apf::Mesh* am,
    osh_t om)
{
  unsigned nf = osh_nfields(om, 0);
  for (unsigned i = 0; i < nf; ++i) {
    std::string name = osh_field(om, 0, i);
    if (name != "coordinates" && name != "param")
      fieldToAPF(am, om, name.c_str());
  }
}

static void coordinatesFromAPF(
    apf::Mesh* am,
    osh_t om)
{
  fieldFromAPF(am, om, am->getCoordinateField());
}

static void parametricFromAPF(
    apf::Mesh* am,
    osh_t om)
{
  osh_new_field(om, 0, "param", 3, OSH_TRANSFER_NOT);
  double* param = osh_get_field(om, 0, "param");
  apf::MeshEntity* v;
  apf::MeshIterator* it = am->begin(0);
  unsigned i = 0;
  while ((v = am->iterate(it))) {
    apf::Vector3 pt;
    am->getParam(v, pt);
    pt.toArray(param + i * 3);
    ++i;
  }
  am->end(it);
}

osh_t fromAPF(apf::Mesh* am)
{
  apf::Numbering* local = apf::numberOverlapNodes(am, "osh_id");
  unsigned dim = (unsigned) am->getDimension();
  osh_t om = osh_new(dim);
  unsigned nverts = am->count(0);
  osh_build_ents(om, 0, nverts);
  coordinatesFromAPF(am, om);
  parametricFromAPF(am, om);
  for (unsigned d = 1; d <= dim; ++d)
    connectivityFromAPF(am, om, local, d);
  for (unsigned d = 0; d <= dim; ++d) {
    classificationFromAPF(am, om, d);
    globalFromAPF(am, om, d);
  }
  fieldsFromAPF(am, om);
  return om;
}

static apf::MeshEntity** verticesToAPF(
    osh_t om,
    apf::Mesh2* am)
{
  apf::MeshEntity** verts = new apf::MeshEntity*[osh_nverts(om)];
  double const* coords = osh_coords(om);
  double const* param = osh_get_field(om, 0, "param");
  unsigned const* class_dim = osh_get_label(om, 0, "class_dim");
  unsigned const* class_id = osh_get_label(om, 0, "class_id");
  for (unsigned i = 0; i < osh_nverts(om); ++i) {
    apf::ModelEntity* ge = am->findModelEntity(class_dim[i], class_id[i]);
    apf::Vector3 x(coords + i * 3);
    apf::Vector3 p(param + i * 3);
    verts[i] = am->createVertex(ge, x, p);
  }
  assert(am->count(0) == osh_nverts(om));
  return verts;
}

static apf::MeshEntity** entitiesToAPF(
    osh_t om,
    apf::Mesh2* am,
    apf::MeshEntity* const* verts,
    unsigned ent_dim)
{
  apf::MeshEntity** ents = new apf::MeshEntity*[osh_count(om, ent_dim)];
  unsigned const* down = osh_down(om, ent_dim, 0);
  unsigned const* class_dim = osh_get_label(om, ent_dim, "class_dim");
  unsigned const* class_id = osh_get_label(om, ent_dim, "class_id");
  apf::Mesh::Type t = apf::Mesh::simplexTypes[ent_dim];
  unsigned nverts_per_ent = ent_dim + 1;
  for (unsigned i = 0; i < osh_count(om, ent_dim); ++i) {
    apf::ModelEntity* ge = am->findModelEntity(class_dim[i], class_id[i]);
    apf::Downward ev;
    for (unsigned j = 0; j < nverts_per_ent; ++j)
      ev[j] = verts[down[i * nverts_per_ent + j]];
    ents[i] = apf::buildElement(am, ge, t, ev);
  }
  return ents;
}

static void ownersToAPF(
    osh_t om,
    apf::Mesh2* am,
    apf::MeshEntity* const* ents,
    unsigned ent_dim)
{
  unsigned const* own_ranks = osh_own_rank(om, ent_dim);
  unsigned const* own_ids = osh_own_id(om, ent_dim);
  /* currently MDS defines ownership between pairs of ranks
     only, which is insufficient to represent ghosting
     ownership of elements.
     therefore we will tag the owners to the entities,
     which can later be used by an apf::Sharing */
  apf::MeshTag* own_tag = am->findTag("owner");
  if (!own_tag)
    own_tag = am->createIntTag("owner", 1);
  for (unsigned i = 0; i < osh_count(om, ent_dim); ++i) {
    int yay_for_signed = own_ranks[i];
    am->setIntTag(ents[i], own_tag, &yay_for_signed);
  }
  PCU_Comm_Begin();
  for (unsigned i = 0; i < osh_count(om, ent_dim); ++i) {
    PCU_COMM_PACK(own_ranks[i], own_ids[i]);
    PCU_COMM_PACK(own_ranks[i], ents[i]);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    unsigned own_id;
    PCU_COMM_UNPACK(own_id);
    apf::MeshEntity* r;
    PCU_COMM_UNPACK(r);
    int from = PCU_Comm_Sender();
    if (from == PCU_Comm_Self()) 
    {
      if ((unsigned)from != own_ranks[own_id])
        printf("fatal error - from %d, own_ranks[own_id] %d\n", from, own_ranks[own_id]);
      assert((unsigned)from == own_ranks[own_id]);
      continue;
    }
    apf::MeshEntity* owner = ents[own_id];
    am->addRemote(owner, from, r);
  }
  PCU_Comm_Begin();
  for (unsigned i = 0; i < osh_count(om, ent_dim); ++i)
    if (own_ranks[i] == (unsigned)(PCU_Comm_Self())) {
      apf::Copies remotes;
      am->getRemotes(ents[i], remotes);
      unsigned ncopies = remotes.size();
      int self = PCU_Comm_Self();
      APF_ITERATE(apf::Copies, remotes, it) {
        PCU_COMM_PACK(it->first, it->second);
        PCU_COMM_PACK(it->first, ncopies);
        PCU_COMM_PACK(it->first, self);
        PCU_COMM_PACK(it->first, ents[i]);
        APF_ITERATE(apf::Copies, remotes, it2)
          if (it2->first != it->first) {
            PCU_COMM_PACK(it->first, it2->first);
            PCU_COMM_PACK(it->first, it2->second);
          }
      }
    }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    apf::MeshEntity* e;
    PCU_COMM_UNPACK(e);
    unsigned ncopies;
    PCU_COMM_UNPACK(ncopies);
    for (unsigned i = 0; i < ncopies; ++i) {
      int p;
      PCU_COMM_UNPACK(p);
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      am->addRemote(e, p, r);
    }
  }
}

void toAPF(osh_t om, apf::Mesh2* am)
{
  apf::MeshEntity** ents[4];
  ents[0] = verticesToAPF(om, am);
  for (unsigned d = 1; d <= osh_dim(om); ++d)
    ents[d] = entitiesToAPF(om, am, ents[0], d);
  for (unsigned d = 0; d <= osh_dim(om); ++d) {
    ownersToAPF(om, am, ents[d], d);
    apf::initResidence(am, d);
  }
  for (unsigned d = 0; d <= osh_dim(om); ++d)
    delete [] ents[d];
  am->acceptChanges();
  fieldsToAPF(am, om);
}

};
