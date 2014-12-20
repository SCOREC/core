#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include <apf.h>
#include <limits.h>
#include "maximalIndependentSet/mis.h"


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
        PCU_Debug_Open();
        apf::Parts res;
        const double t1 = MPI_Wtime();
        bool isMis = isInMIS(m);
        double elapsedTime = MPI_Wtime() - t1;
        PCU_Max_Doubles(&elapsedTime, 1);
        if( !PCU_Comm_Self() )
          fprintf(stdout,"mis completed in %f (seconds)\n", elapsedTime);
        
        int side = 0;
        PCU_Comm_Begin();
        if(getSmallSide(s, 0.5*getAvgSides(s), side) && isMis) {
          PCU_Comm_Pack(side,NULL,0);
          getOtherRes(m, s, side, res);
        }

        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          side = PCU_Comm_Sender();
          getOtherRes(m, s, side, res);
        }
        PCU_Comm_Begin();
        PCU_Debug_Print("res ");
        APF_ITERATE(apf::Parts, res, r) {
          PCU_Comm_Pack(*r, NULL, 0);
          PCU_Debug_Print(" %d ", *r);
        }
        PCU_Debug_Print("\n");
        PCU_Comm_Send();
        while (PCU_Comm_Listen())
          setTarget(PCU_Comm_Sender(), s, w, alpha);
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
      bool isInMIS(apf::Mesh* m) {
        misLuby::partInfo part;
        part.id = PCU_Comm_Self();
        part.net.push_back(PCU_Comm_Self());
  
        apf::Parts neighbors;
        apf::MeshEntity* vtx;
        apf::MeshIterator* vtxs = m->begin(0);
        while ((vtx = m->iterate(vtxs))) {
          apf::Parts residence;
          m->getResidence(vtx,residence);
          apf::unite(neighbors,residence);
        }
        m->end(vtxs);
        neighbors.erase(m->getId());

        for (apf::Parts::iterator itr = neighbors.begin();
             itr!=neighbors.end();itr++) {
          part.adjPartIds.push_back(*itr);
          part.net.push_back(*itr);
        }
        const int t = static_cast<int>(MPI_Wtime());
        int randNumSeed = t+PCU_Comm_Self()+1;
        mis_init(randNumSeed,true);
        bool isIn =mis(part, false, true);
        return isIn;
      }
  };
  Targets* makeShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha) {
    return new ShapeTargets(m,s,w,alpha);
  }
} //end namespace
