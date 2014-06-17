#include "parma_base.h"
#include "parma_associative.h"
#include <PCU.h>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <string>
#include <unistd.h>

namespace parma {
  class Sides : public Associative<int> {
    public:
      Sides(apf::Mesh* m) {
        totalSides = 0;
        init(m);
      }
      int total() {
        return totalSides;
      }
    private:
      int totalSides;
      void init(apf::Mesh* m) {
        apf::MeshEntity* s;
        apf::MeshIterator* it = m->begin(m->getDimension()-1);
        totalSides = 0;
        while ((s = m->iterate(it)))
          if (m->countUpward(s)==1 && m->isShared(s)) {
            const int peerId = apf::getOtherCopy(m,s).first;
            set(peerId, get(peerId)+1);
            ++totalSides;
          }
        m->end(it);
      }
  };

  class ElmTargets : public Targets {
    public:
      ElmTargets(Sides* s, Weights* w, double alpha) {
        init(s, w, alpha);
      }
      double total() {
        return totW;
      }
    private:
      double totW;
      void init(Sides* s, Weights* w, double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double selfW = w->self();
          const double peerW = w->get(peer);
          if ( selfW > peerW ) {
            const double difference = selfW - peerW;
            double sideFraction = side->second;
            sideFraction /= s->total();
            double peerW = difference * sideFraction * alpha;
            set(peer, peerW);
            totW+=peerW;
          }
        }
        s->end();
      }
  };

  Balancer::Balancer(apf::Mesh* mIn, apf::MeshTag* wIn, double alphaIn) 
    : m(mIn), w(wIn), alpha(alphaIn) {
      sides = new Sides(m);
  }

  void Balancer::setSelector(Selector* s) { selects = s; }
  void Balancer::setWeights(Weights* w) { weights = w };
  void Balancer::setTargets(Targets* t) { targets = t };

  Balancer::~Balancer() {
    delete sides;
    delete weights;
    delete targets;
    delete selects;
  }

  bool Balancer::run(double maxImb, int verbosityIn) {
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
