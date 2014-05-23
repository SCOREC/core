#include "parma_rib.h"
#include <apfPartition.h>

namespace parma {

static apf::Migration* splitMesh(apf::Mesh* m, apf::MeshTag* weights, int depth)
{
  int dim = m->getDimension();
  apf::DynamicArray<Body> arr(m->count(dim));
  apf::DynamicArray<apf::MeshEntity*> elems(m->count(dim));
  apf::MeshEntity* e;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it))) {
    arr[i].point = getLinearCentroid(m, e);
    if (weights)
      m->getDoubleTag(e, weights, &(arr[i].mass));
    else
      arr[i].mass = 1;
    elems[i] = e;
    ++i;
  }
  assert(i == m->count(dim));
  m->end(it);
  Bodies all(&arr[0], arr.getSize());
  int n = 1 << depth;
  apf::DynamicArray<Bodies> out(n);
  recursivelyBisect(&all, depth, &out[0]);
  apf::Migration* plan = new apf::Migration(m);
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < out[i].n; ++j) {
      size_t k = out[i].body[j] - &arr[0];
      plan->send(elems[k], i);
    }
  }
  all.destroy();
  return plan;
}

class RibSplitter : public apf::Splitter
{
  public:
    RibSplitter(apf::Mesh* m)
    {
      mesh = m;
    }
    virtual ~RibSplitter() {}
    virtual apf::Migration* split(apf::MeshTag* weights, double tolerance,
        int multiple)
    {
      int depth;
      for (depth = 0; (1 << depth) < multiple; ++depth);
      assert((1 << depth) == multiple);
      return splitMesh(mesh, weights, depth);
    }
  private:
    apf::Mesh* mesh;
};

}

apf::Splitter* Parma_MakeRibSplitter(apf::Mesh* m)
{
  return new parma::RibSplitter(m);
}

