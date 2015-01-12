#include <PCU.h>
#include "maBalance.h"
#include "maAdapt.h"
#include <parma.h>
#include <apfZoltan.h>

namespace ma {

static double clamp(double x, double max, double min)
{
  if (x > max) return max;
  if (x < min) return min;
  return x;
}

static double getSizeWeight(Adapt* a, Entity* e, int type)
{
/* until we have a size field that matches the
   anisotropy of the layer itself, we have to
   hack around to prevent prisms from getting
   very small weights due to their very small thickness.
   the current hack will be to measure their
   triangular bases instead.
   pyramids will get very small weights, but
   they are supposed to be fairly sparse so
   that should not ruin things. */
  if (type == PRISM) {
    Entity* f[5];
    a->mesh->getDownward(e, 2, f);
    return a->sizeField->getWeight(f[0]);
  }
  return a->sizeField->getWeight(e);
}

static double clampForIterations(Adapt* a, double weight)
{
  Mesh* m = a->mesh;
  int dimension = m->getDimension();
  double max = pow(2.0, dimension*(a->refinesLeft));
/* coarsening performance is more empirical: 3x decrease in tet
   count when uniformly refining a 58k element cube, 4x decrease
   on some 2D meshes which were more structured. */
  double min = pow(4.0, -(a->coarsensLeft));
  return clamp(weight, max, min);
}

static double clampForLayerPermissions(Adapt* a, int type, double weight)
{
  if (apf::isSimplex(type))
    return weight;
  if ( ! a->input->shouldRefineLayer)
    weight = std::max(1.0, weight);
  if ( ! a->input->shouldCoarsenLayer)
    weight = std::min(1.0, weight);
  return weight;
}

static double accountForTets(Adapt* a, int type, double weight)
{
  if (a->input->shouldTurnLayerToTets) {
    if (type == PRISM)
      return weight * 3;
    if (type == PYRAMID)
      return weight * 2;
  }
  return weight;
}

double getElementWeight(Adapt* a, Entity* e)
{
  int type = a->mesh->getType(e);
  double weight = getSizeWeight(a, e, type);
  weight = clampForIterations(a, weight);
  weight = clampForLayerPermissions(a, type, weight);
  return accountForTets(a, type, weight);
}

Tag* getElementWeights(Adapt* a)
{
  Mesh* m = a->mesh;
  Tag* weights = m->createDoubleTag("ma_weight",1);
  Entity* e;
  Iterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
  {
    double weight = getElementWeight(a,e);
    m->setDoubleTag(e,weights,&weight);
  }
  m->end(it);
  return weights;
}

static void runBalancer(Adapt* a, apf::Balancer* b)
{
  Mesh* m = a->mesh;
  Input* in = a->input;
  Tag* weights = getElementWeights(a);
  b->balance(weights,in->maximumImbalance);
  delete b;
  removeTagFromDimension(m,weights,m->getDimension());
  m->destroyTag(weights);
}

void runZoltan(Adapt* a)
{
  runBalancer(a, apf::makeZoltanBalancer(
        a->mesh, apf::GRAPH, apf::REPARTITION,
        /* debug = */ false));
}

void runParma(Adapt* a)
{
  runBalancer(a, Parma_MakeCentroidDiffuser(a->mesh));
}

void printEntityImbalance(Mesh* m)
{
  double imbalance[4];
  Parma_GetEntImbalance(m,&imbalance);
  double p = (imbalance[m->getDimension()]-1)*100;
  print("element imbalance %.0f%% of average",p);
}

void preBalance(Adapt* a)
{
  if (PCU_Comm_Peers()==1)
    return;
  Input* in = a->input;
  if (in->shouldRunPreZoltan)
    runZoltan(a);
  if (in->shouldRunPreParma)
    runParma(a);
}

void midBalance(Adapt* a)
{
  if (PCU_Comm_Peers()==1)
    return;
  Input* in = a->input;
  if (in->shouldRunMidZoltan)
    runZoltan(a);
  if (in->shouldRunMidParma)
    runParma(a);
}

void postBalance(Adapt* a)
{
  if (PCU_Comm_Peers()==1)
    return;
  Input* in = a->input;
  if (in->shouldRunPostZoltan)
    runZoltan(a);
  if (in->shouldRunPostParma)
    runParma(a);
  printEntityImbalance(a->mesh);
}

}
