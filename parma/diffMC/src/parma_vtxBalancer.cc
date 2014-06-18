namespace parma {
  class VtxBalancer : public apf::Balancer {
    public:
      VtxBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
        }
      bool runStep(apf::MeshTag* weights, double tolerance) {
        parma::Balancer parmaVtx(mesh, weights, layers, bridge, factor);
        parmaVtx.setSides(makeElmBdrySides(m));
        parmaVtx.setSelector(new VtxSelector(m, w));
        parmaVtx.setWeights(new EntWeights(m, w, SIDES, m->getDimension()));
        parmaVtx.setTargets(new ElmTargets(CAKE, CAKES, CATs));
        return parmaVtx.run(tolerance);
      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        double t0 = MPI_Wtime();
        while (runStep(weights,tolerance));
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("vertices balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
}; //end parma namespace

apf::Balancer* Parma_MakeVtxBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::VtxBalancer(m, stepFactor, verbosity);
}
