#include "apfZoltan.h"
#include "apfZoltanBridge.h"
#include <apfPartition.h>
#include <PCU.h>

namespace apf {

class ZoltanSplitter : public Splitter
{
  public:
    ZoltanSplitter(Mesh* m, int method, int approach, bool sync):
      bridge(m, true, method, approach)
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
    ZoltanBridge bridge;
};

class ZoltanBalancer : public Balancer
{
  public:
    ZoltanBalancer(Mesh* m, int method, int approach):
      bridge(m, false, method, approach)
    {}
    virtual ~ZoltanBalancer() {}
    virtual void balance(MeshTag* weights, double tolerance)
    {
      Migration* plan = bridge.run(weights, tolerance, 1);
      bridge.mesh->migrate(plan);
    }
  private:
    ZoltanBridge bridge;
};

Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach, bool sync)
{
  return new ZoltanSplitter(mesh, method, approach, sync);
}

Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach)
{
  return new ZoltanBalancer(mesh, method, approach);
}

}
