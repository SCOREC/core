#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include <apf.h>
#include <limits.h>
namespace parma {
  class ShapeTargets : public Targets {
    public:
      ShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha) {
        init(m,s,w,alpha);
      }
      double total() {
        return totW;
      }
    private:
      ShapeTargets();
      double totW;
      void init(apf::Mesh* m, Sides* s, Weights* w, double alpha) {
        apf::Parts res;
        int side = 0;
        if( getSmallSide(s, 0.5*getAvgSides(s), side) )
          getOtherRes(m, s, side, res);
        PCU_Comm_Begin();
        PCU_Debug_Print("res ");
        APF_ITERATE(apf::Parts, res, r) {
          PCU_Comm_Pack(*r, NULL, 0);
          PCU_Debug_Print(" %d ", *r);
        }
        PCU_Debug_Print("\n");
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) 
          if( ! res.count(PCU_Comm_Sender()) )
            setTarget(PCU_Comm_Sender(), s, w, alpha);
          else 
            PCU_Debug_Print("target conflict res has %d\n", PCU_Comm_Sender());
      }
      void getOtherRes(apf::Mesh* m, Sides*, int peer, apf::Parts& res) {
        const int self = PCU_Comm_Self();
        apf::MeshEntity* e;
        apf::MeshIterator* itr = m->begin(m->getDimension()-2);
        while( (e = m->iterate(itr)) ) {
          apf::Parts eRes;
          m->getResidence(e, eRes);
          if( eRes.count(self) && eRes.count(peer) ) 
            APF_ITERATE(apf::Parts, eRes, r)
              res.insert(*r);
        }
        m->end(itr);
        res.erase(peer);
        res.erase(self);
      }
      double getAvgSides(Sides* s) {
        double tot = s->total();
        PCU_Add_Doubles(&tot, 1);
        int cnt = s->size();
        PCU_Add_Ints(&cnt, 1);
        return tot/cnt;
      }
      bool getSmallSide(Sides* s, double small, int& peer) {
        int minSides = INT_MAX;
        peer = -1;
        s->begin();
        const Sides::Item* side;
        while( (side = s->iterate()) ) 
          if( side->second < small && side->second < minSides ) {
            peer = side->first;
            minSides = side->second;
          }
        s->end();
        return (peer != -1);
      }
      void setTarget(const int peer, Sides* s, Weights* w, double alpha) {
        assert(s->has(peer));
        double sideFactor = (double)s->get(peer) / s->total();
        double scaledW = alpha * sideFactor * w->self();
        set(peer, scaledW);
        totW+=scaledW;
      }
  };
  Targets* makeShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha) {
    return new ShapeTargets(m,s,w,alpha);
  }
} //end namespace
