#include "PCU.h"
#include "parma_dcpart.h"
#include "parma_commons.h"
#include <maximalIndependentSet/mis.h>

#define TO_UINT(a) static_cast<unsigned>(a)
#define TO_INT(a) static_cast<int>(a)

typedef std::map<unsigned, unsigned> muu;

namespace {
  bool isInMis(muu& mt) {
    unsigned seed = TO_UINT(PCU_Comm_Self()+1);
    mis_init(seed);
    misLuby::partInfo part;
    part.id = PCU_Comm_Self();
    std::set<int> targets;
    APF_ITERATE(muu, mt, mtItr) {
      int peer = TO_INT(mtItr->second);
      if( !targets.count(peer) ) {
        part.adjPartIds.push_back(peer);
        part.net.push_back(peer);
        targets.insert(peer);
      }
    }
    part.net.push_back(part.id);
    return mis(part,false,true);
  }
}

class dcPartFixer::PartFixer : public dcPart {
  public:
    PartFixer(apf::Mesh* mesh, unsigned verbose=0) 
      : dcPart(mesh,verbose), m(mesh), vb(verbose) 
    {
      fix();
    }

  private:
    apf::Mesh* m;
    unsigned vb;

    int totNumDc() {
      int ndc = TO_INT(numDisconnectedComps());
      PCU_Add_Ints(&ndc, 1);
      return ndc;
    }

    void setupPlan(muu& dcCompTgts, apf::Migration* plan) {
      apf::MeshEntity* e;
      apf::MeshIterator* itr = m->begin(m->getDimension());
      while( (e = m->iterate(itr)) ) {
        if( isIsolated(e) ) continue;
        unsigned id = compId(e);
        if ( dcCompTgts.count(id) )
          plan->send(e, TO_INT(dcCompTgts[id]));
      }
      m->end(itr);
    }

    /**
     * @brief remove the disconnected set(s) of elements from the part
     * @remark migrate the disconnected set(s) of elements into the adjacent part
     *         that shares the most faces with the disconnected set of elements
     *         requires that the sets of elements forming disconnected components
     *         are tagged
     */
    void fix() {
      double t1 = PCU_Time();
      int loop = 0;
      int ndc = 0;
      while( (ndc = totNumDc()) && loop++ < 50 ) {
        double t2 = PCU_Time();
        muu dcCompTgts;

        unsigned maxSz = 0;
        for(unsigned i=0; i<getNumComps(); i++)
          if( getCompSize(i) > maxSz )
            maxSz = getCompSize(i);

        for(unsigned i=0; i<getNumComps(); i++)
          if( getCompSize(i) != maxSz )
            dcCompTgts[i] = getCompPeer(i);
        assert( dcCompTgts.size() == getNumComps()-1 );
        apf::Migration* plan = new apf::Migration(m);
        if ( isInMis(dcCompTgts) )
          setupPlan(dcCompTgts, plan);

        reset();
        double t3 = PCU_Time();
        m->migrate(plan);
        if( ! PCU_Comm_Self() && vb)
          parmaCommons::status(
              "loop %d components %d seconds <fix migrate> %.3f %.3f\n",
              loop, ndc, t3-t2, PCU_Time()-t3);
      }
      parmaCommons::printElapsedTime(__func__, PCU_Time() - t1);
    }
};

dcPartFixer::dcPartFixer(apf::Mesh* mesh, unsigned verbose) 
: pf( new PartFixer(mesh,verbose) ) {}

dcPartFixer::~dcPartFixer() {
delete pf;
}
