#include <PCU.h>
#include <apfPartition.h>

namespace {
  void printTiming(const char* type, int steps, double tol, double time) {
    if (!PCU_Comm_Self())
      printf("%s balanced in %d steps to %f in %f seconds\n",
          type, steps, tol, time);
  }
}

namespace parma {
  class Balancer : public apf::Balancer {
    public:
      Balancer(apf::Mesh* m, double f, int v, const char* n)
        : mesh(m), factor(f), verbose(v), name(n) {
        maxStep = 50;
      }
      virtual bool runStep(apf::MeshTag* wtag, double tolerance)=0;
      virtual void balance(apf::MeshTag* wtag, double tolerance) {
        int step = 0;
        double t0 = PCU_Time();
        while (runStep(wtag,tolerance) && step++ < maxStep);
        printTiming(name, step, tolerance, PCU_Time()-t0);
      }
      apf::Mesh* mesh;
      double factor;
      int verbose;
      const char* name;
      int maxStep;
  };
}
