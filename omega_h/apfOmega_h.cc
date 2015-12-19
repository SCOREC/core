#include "apfOmega_h.h"

#include <cstdlib>

#include <apfMesh2.h>
#include <apfNumbering.h>

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

osh_t fromAPF(apf::Mesh* am)
{
  apf::Numbering* local = apf::numberOverlapNodes(am, "osh_id");
  unsigned dim = (unsigned) am->getDimension();
  osh_t om = osh_new(dim);
  unsigned nverts = am->count(0);
  osh_build_ents(om, 0, nverts);
  connectivityFromAPF(am, om, local, dim);
  apf::GlobalNumbering* glob_n = apf::makeGlobal(
      apf::numberOwnedNodes(am, "osh_global"));
  osh_new_field(om, "coordinates", 3);
  double* coords = osh_get_field(om, "coordinates");
  osh_new_label(om, "class_dim", 1);
  unsigned* class_dim = osh_get_label(om, "class_dim");
  osh_new_label(om, "class_id", 1);
  unsigned* class_id = osh_get_label(om, "class_id");
  unsigned long* global = (unsigned long*) malloc(
      sizeof(unsigned long) * nverts);
  {
    apf::MeshEntity* v;
    apf::MeshIterator* it = am->begin(0);
    unsigned i = 0;
    while ((v = am->iterate(it))) {
      apf::Vector3 pt;
      am->getPoint(v, 0, pt);
      pt.toArray(coords + i * 3);
      apf::ModelEntity* ge = am->toModel(v);
      class_dim[i] = am->getModelType(ge);
      class_id[i] = am->getModelTag(ge);
      global[i] = apf::getNumber(glob_n, apf::Node(v, 0));
      ++i;
    }
    am->end(it);
  }
  apf::destroyGlobalNumbering(glob_n);
  return om;
}

apf::Mesh2* toAPF(osh_t om)
{
  (void) om;
  /* todo ? todo. */
  return 0;
}

};
