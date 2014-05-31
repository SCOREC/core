#include "maBalance.h"
#include "maAdapt.h"
#include <parma.h>
#include <apfZoltan.h>
#include <PCU.h>

namespace ma {

double clamp(double x, double max, double min)
{
  if (x > max) return max;
  if (x < min) return min;
  return x;
}

double getSimplexWeight(Adapt* a, Entity* e)
{
  Mesh* m = a->mesh;
  SizeField* sf = a->sizeField;
  int dimension = m->getDimension();
  double max = pow(2.0, dimension*(a->refinesLeft));
/* coarsening performance is more empirical: 3x decrease in tet
   count when uniformly refining a 58k element cube, 4x decrease
   on some 2D meshes which were more structured. we will be conservative
   towards more elements
   since light parts are better than heavy parts */
  double min = pow(3.0, -(a->coarsensLeft));
/* get the measurement in metric space - best estimate of how many tets
   this tet would turn into after infinite iterations */
  double weight = sf->getWeight(e);
/* and clamp it based on how much we could actually refine or coarsen
   this tet in the remaining iterations */
  return clamp(weight,max,min);
}

double getLayerWeight(Adapt* a, Entity* e)
{
  if (a->input->shouldTurnLayerToTets)
  {
    Mesh* m = a->mesh;
    int type = m->getType(e);
    if (type==PRISM)
      return 3.0;
    else
    { assert(type==PYRAMID);
      return 2.0;
    }
  }
  return 1.0;
}

double getElementWeight(Adapt* a, Entity* e)
{
  if (apf::isSimplex(a->mesh->getType(e)))
    return getSimplexWeight(a,e);
  else
    return getLayerWeight(a,e);
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
        a->mesh, apf::GRAPH, apf::REPARTITION));
}

void runDiffusion(Adapt* a)
{
  runBalancer(a, Parma_MakeCentroidDiffuser(a->mesh));
}

void runParma(Adapt* a)
{
  double t0 = MPI_Wtime();
  Mesh* m = a->mesh;
  Input* in = a->input;
  Tag* weights = getElementWeights(a);
  double entityImbalance[4];
  Parma_GetWeightedEntImbalance(m,weights,&entityImbalance);
  double& elementImbalance = entityImbalance[m->getDimension()];
  print("element imbalance before parma %f",elementImbalance);
  int priorities[4] = {0,0,0,1};
  Parma_RunWeightedPtnImprovement(
      m,
      weights,
      &priorities,
      in->maximumImbalance,
      0,
      in->diffuseIterations);
  removeTagFromDimension(m,weights,m->getDimension());
  m->destroyTag(weights);
  double t1 = MPI_Wtime();
  print("parma run took %f seconds",t1-t0);
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
  else if (in->shouldRunPreParma)
    runParma(a);
  else if (in->shouldRunPreDiffusion)
    runDiffusion(a);
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
  else if (in->shouldRunMidDiffusion)
    runDiffusion(a);
}

void postBalance(Adapt* a)
{
  if (PCU_Comm_Peers()==1)
    return;
  Input* in = a->input;
  if (in->shouldRunPostZoltan)
    runZoltan(a);
  if (in->shouldRunPostDiffusion)
    runDiffusion(a);
  else if (in->shouldRunPostParma)
    runParma(a);
  printEntityImbalance(a->mesh);
}

}
