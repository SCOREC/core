#include "parma_base.h"
#include <apf_partition.h>

namespace parma {
  class HpsTargets : public Targets {
    public:
      HpsTargets(Sides* s, Weights* w) 
        : Targets(s,w,0.1)
      {
        init(s, w);
      }
      double total() {
        return totW;
      }
    private:
      HpsTargets();
      double totW;
      void init(Sides* s, Weights* w) {
      }
  };

  class HpsSelector : public Selector {
    public:
      HpsSelector(apf::Mesh* m, apf::MeshTag* w)
        : Selector(m, w) {}
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        //determine where elements need to migrated (merges and splits) 
        // and insert them into the plan
        return plan;
      }
    private:
      HpsSelector();
      apf::Mesh* mesh;
      apf::MeshTag* vtag;
      apf::MeshTag* wtag;
  };

  class HpsBalancer : public apf::Balancer {
    public:
      HpsBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
        }
      bool run(apf::MeshTag* weights, double tolerance) {
        parma::Balancer hps(mesh, weights, layers, bridge, factor);
        Sides* sides = makeElmBdrySides(m);
        parmaVtx.setSides(sides);
        Weights* weights = makeEntWeights(m, w, sides, m->getDimension());
        hps.setWeights(weights);
        Targets* targets = makeTargets(sides, weights, factor);
        hps.setTargets(targets);
        hps.setSelector(new HpsSelector(m, w));
        return hps.run(tolerance);
      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        double t0 = MPI_Wtime();
        run(weights,tolerance);
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("elements balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
}; //end parma namespace

apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::HpsBalancer(m, stepFactor, verbosity);
}
