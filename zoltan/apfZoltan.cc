/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfZoltan.h"
#include "apfZoltanMesh.h"
#include <apfPartition.h>
#include <lionPrint.h>

namespace apf {

class ZoltanSplitter : public Splitter
{
  public:
    ZoltanSplitter(Mesh* m, int method, int approach,
        bool debug, bool sync, bool local=true):
      bridge(m, local, method, approach, debug)
    {
      isSynchronous = sync;
    }
    virtual ~ZoltanSplitter() {}
    virtual Migration* split(MeshTag* weights, double tolerance, int multiple)
    {
      double t0 = pcu::Time();
      Migration* plan = bridge.run(weights, tolerance, multiple);
      if (isSynchronous) {
        for (int i = 0; i < plan->count(); ++i) {
          MeshEntity* e = plan->get(i);
          int p = plan->sending(e);
          p += bridge.mesh->getPCU()->Self() * multiple;
          plan->send(e, p);
        }
      }
      double t1 = pcu::Time();
      if (!bridge.mesh->getPCU()->Self())
        lion_oprint(1, "planned Zoltan split factor %d to target"
            " imbalance %f in %f seconds\n", multiple, tolerance, t1 - t0);
      return plan;
    }
  private:
    bool isSynchronous;
    ZoltanMesh bridge;
};

class ZoltanBalancer : public Balancer
{
  public:
    ZoltanBalancer(Mesh* m, int method, int approach, bool debug):
      bridge(m, false, method, approach, debug)
    {}
    virtual ~ZoltanBalancer() {}
    virtual void balance(MeshTag* weights, double tolerance)
    {
      double t0 = pcu::Time();
      Migration* plan = bridge.run(weights, tolerance, 1);
      if (!bridge.mesh->getPCU()->Self())
        lion_oprint(1, "planned Zoltan balance to target "
            "imbalance %f in %f seconds\n",
            tolerance, pcu::Time() - t0);
      bridge.mesh->migrate(plan);
      double t1 = pcu::Time();
      if (!bridge.mesh->getPCU()->Self())
        lion_oprint(1,"Zoltan balanced to %f in %f seconds\n",
            tolerance, t1-t0);
    }
  private:
    ZoltanMesh bridge;
};

Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach,
    bool debug, bool sync)
{
  return new ZoltanSplitter(mesh, method, approach, debug, sync);
}

Splitter* makeZoltanGlobalSplitter(Mesh* mesh, int method, int approach,
    bool debug) {
  bool sync = false; //don't alter the plan
  bool local = false;
  return new ZoltanSplitter(mesh, method, approach, debug, sync, local);
}

Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach,
    bool debug)
{
  return new ZoltanBalancer(mesh, method, approach, debug);
}

}
