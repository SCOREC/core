#include "parma_dcpart.h"
#include "parma_commons.h"
#include "parma_convert.h"
#include <maximalIndependentSet/mis.h>
#include <pcu_util.h>

typedef std::map<unsigned, unsigned> muu;

namespace {
  bool isInMis(muu& mt, pcu::PCU *PCUObj) {
    unsigned seed = TO_UINT(PCUObj->Self()+1);
    mis_init(seed);
    misLuby::partInfo part;
    part.id = PCUObj->Self();
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
    return mis(part,PCUObj,false,true);
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
      return m->getPCU()->Add<int>(ndc);
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
      double t1 = pcu::Time();
      int loop = 0;
      int ndc = 0;
      while( (ndc = totNumDc()) && loop++ < 50 ) {
        double t2 = pcu::Time();
        muu dcCompTgts;

        unsigned maxSz = 0;
        for(unsigned i=0; i<getNumComps(); i++)
          if( getCompSize(i) > maxSz )
            maxSz = getCompSize(i);

        for(unsigned i=0; i<getNumComps(); i++)
          if( getCompSize(i) != maxSz )
            dcCompTgts[i] = getCompPeer(i);
        PCU_ALWAYS_ASSERT( dcCompTgts.size() == getNumComps()-1 );
        apf::Migration* plan = new apf::Migration(m);
        if ( isInMis(dcCompTgts, m->getPCU()) )
          setupPlan(dcCompTgts, plan);

        reset();
        double t3 = pcu::Time();
        m->migrate(plan);
        if( ! m->getPCU()->Self() && vb)
          parmaCommons::status(
              "loop %d components %d seconds <fix migrate> %.3f %.3f\n",
              loop, ndc, t3-t2, pcu::Time()-t3);
      }
      parmaCommons::printElapsedTime(__func__, pcu::Time() - t1, m->getPCU());
    }
};

dcPartFixer::dcPartFixer(apf::Mesh* mesh, unsigned verbose) 
: pf( new PartFixer(mesh,verbose) ) {}

dcPartFixer::~dcPartFixer() {
delete pf;
}
