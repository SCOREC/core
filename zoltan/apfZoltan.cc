#include "apfZoltan.h"
#include "apfZoltanMesh.h"
#include <apfPartition.h>
#include <PCU.h>

namespace apf {

class ZoltanSplitter : public Splitter
{
  public:
    ZoltanSplitter(Mesh* m, int method, int approach, bool debug, bool sync):
      bridge(m, true, method, approach, debug)
    {
      isSynchronous = sync;
    }
    virtual ~ZoltanSplitter() {}
    virtual Migration* split(MeshTag* weights, double tolerance, int multiple)
    {
      Migration* plan = bridge.run(weights, tolerance, multiple);
      if (isSynchronous)
        for (int i = 0; i < plan->count(); ++i) {
          MeshEntity* e = plan->get(i);
          int p = plan->sending(e);
          p += PCU_Proc_Self() * multiple;
          plan->send(e, p);
        }
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
      Migration* plan = bridge.run(weights, tolerance, 1);
      bridge.mesh->migrate(plan);
    }
  private:
    ZoltanMesh bridge;
};

Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach,
    bool debug, bool sync)
{
  return new ZoltanSplitter(mesh, method, approach, debug, sync);
}

Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach,
    bool debug)
{
  return new ZoltanBalancer(mesh, method, approach, debug);
}

}
