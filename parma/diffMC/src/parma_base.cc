  class Balancer {
    public:
      Balancer(apf::Mesh* mIn, apf::MeshTag* wIn, 
          int layersIn, int bridgeIn, double alphaIn) 
        : m(mIn), w(wIn), layers(layersIn), bridge(bridgeIn), alpha(alphaIn)
      {
      }

      ~Balancer();
      bool run(double maxImb);
    private:
      Balancer();
      apf::Mesh* m;
      apf::MeshTag* w;
      double alpha;
      int verbose;
      double imbalance();
      Sides* sides;
      Weights* weights;
      Targets* targets;
      Selector* selects;
  };

  Balancer::~Balancer() {
    delete sides;
    delete weights;
    delete targets;
    delete selects;
  }

  bool Balancer::run(double maxImb) {
    const double imb = imbalance();
    if ( 0 == PCU_Comm_Self() )
      fprintf(stdout, "imbalance %.3f\n", imb);
    if ( imb < maxImb ) 
      return false;
    apf::Migration* plan = selects->run();
    m->migrate(plan);
    return true;
  }

  double Balancer::imbalance() { 
    double maxWeight = 0, totalWeight = 0;
    maxWeight = totalWeight = weights->self();
    PCU_Add_Doubles(&totalWeight,1);
    PCU_Max_Doubles(&maxWeight, 1);
    double averageWeight = totalWeight / PCU_Comm_Peers();
    return maxWeight / averageWeight;
  }
}

