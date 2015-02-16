#include <PCU.h>
#include <parma.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include <apf.h>
#include <limits.h>
#include "maximalIndependentSet/mis.h"



namespace parma {
  class ShapeTargets : public Targets {
    public:
      ShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha,
		   double avgSideMult, bool isInMIS) {
        init(m,s,w,alpha,avgSideMult,isInMIS);
      }
      double total() {
        return totW;
      }
    private:
      ShapeTargets();
      double totW;
      void init(apf::Mesh* m, Sides* s, Weights* w, double alpha,
		double avgSideMult, bool isInMIS) {
        PCU_Debug_Open();
	double avgSide=getAvgSides(s);
        int side = -1;
	int side_length=INT_MAX;
        PCU_Comm_Begin();
        if(isInMIS && getSmallSide(s, avgSideMult*avgSide, side,side_length)) {
          PCU_COMM_PACK(side,side_length);
	  PCU_Debug_Print("small side with %d\n",side);
        }
	side=-1;
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
	  int temp_length;
	  PCU_COMM_UNPACK(temp_length);
	  if (temp_length<side_length) {
            side = PCU_Comm_Sender();
	    side_length=temp_length;
	    PCU_Debug_Print("recv small side from %d with length %d\n",
                side,side_length);
	  }
        }

        apf::Parts res;
	if (side!=-1)
	  getOtherRes(m,s,side,res);

        PCU_Comm_Begin();
        PCU_Debug_Print("res ");
        APF_ITERATE(apf::Parts, res, r) {
          PCU_COMM_PACK(*r, side);
	  PCU_COMM_PACK(*r, side_length);
          PCU_Debug_Print(" %d ", *r);
        }
        PCU_Debug_Print("\n");
	PCU_Comm_Send();
        int min=INT_MAX; std::pair<int,int> minSide;
	while (PCU_Comm_Listen()) {
	  int temp_side;
	  PCU_COMM_UNPACK(temp_side);
	  int temp_length;
	  PCU_COMM_UNPACK(temp_length);
          if( temp_length < min ) { 
	    min=temp_length;
	    minSide = std::make_pair(PCU_Comm_Sender(),temp_side); 
	  }
	  
	}
	if (min!=INT_MAX) {
	  setTarget(minSide.first, s, w, alpha,m);
	  setTarget(minSide.second, s, w, alpha,m);
	}
      }
      void getOtherRes(apf::Mesh* m, Sides*, int peer, apf::Parts& res) {
        const int self = PCU_Comm_Self();
        apf::MeshEntity* e;
        apf::MeshIterator* itr = m->begin(0);
        while( (e = m->iterate(itr)) ) {
          apf::Parts eRes;
          m->getResidence(e, eRes);
          if( eRes.count(self) && eRes.count(peer) ) {
            //PCU_Debug_Print("other res: ");
            APF_ITERATE(apf::Parts, eRes, r) {
	      //PCU_Debug_Print(" %d ",*r);
              res.insert(*r);
	    }
            //PCU_Debug_Print("\n");
          }
        }
        m->end(itr);
        res.erase(peer);
        res.erase(self);
      }
      double getAvgSides(Sides* s) {
        double tot = s->total();
        PCU_Add_Doubles(&tot, 1);
        int cnt = static_cast<int>(s->size());
        PCU_Add_Ints(&cnt, 1);
        return tot/cnt;
      }
    bool getSmallSide(Sides* s, double small, int& peer, int& minSides) {
        minSides = INT_MAX;
        peer = -1;
        s->begin();
        const Sides::Item* side;
        while( (side = s->iterate()) ) 
          if( side->second > 0 &&side->second < small && side->second < minSides ) {
            peer = side->first;
            minSides = side->second;
          }
        s->end();
        return (peer != -1);
      }
    void setTarget(const int peer, Sides* s, Weights* w, double alpha,apf::Mesh* m) {
      PCU_Debug_Print("Sending to %d\n",peer);
      if (!s->has(peer)){
	fprintf(stdout,"Failure by %d with %d\n",PCU_Comm_Self(),peer);
      }
      assert(s->has(peer));
      const double totSides = static_cast<double>(s->total());
      const double sideFactor = s->get(peer) / totSides;
      double scaledW = alpha * sideFactor * w->self();
      set(peer, scaledW);
      totW+=scaledW;
    }
  };
  Targets* makeShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha,
			    double avgSideMult, bool isInMIS) {
    return new ShapeTargets(m,s,w,alpha,avgSideMult,isInMIS);
  }
} //end namespace
