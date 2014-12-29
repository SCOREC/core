#include <PCU.h>
#include "parma_selector.h"
#include "parma_sides.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <vector>
#include <limits.h>

static int WeldSelectorCalls = 0;

namespace {
  typedef std::set<apf::MeshEntity*> SetEnt;

  typedef std::set<int> SetInt;
  void print(const char* key, SetInt& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(SetInt, d, sItr)
      ss << *sItr << " ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::map<int,int> Mii;

  typedef std::pair<int, apf::MeshEntity*> Copy;
  typedef std::multimap<int, Copy > MMiCopy;
  void print(const char* key, MMiCopy& d) {
    Mii cnts;
    APF_ITERATE(MMiCopy, d, ditr)
      cnts[ditr->first]++;
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(Mii, cnts, citr) 
      ss << "(" << citr->first << "," << citr->second << ") ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::multimap<int,apf::MeshEntity*> MMCopy;
  void print(const char* key, MMCopy& d) {
    Mii cnts;
    APF_ITERATE(MMCopy, d, c)
      cnts[c->first]++;
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(Mii, cnts, c) 
      ss << "(" << c->first << "," << c->second << ") ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  void insert(MMiCopy& d, int peer, int rmt, apf::MeshEntity* e) {
    Copy c = std::make_pair(rmt,e);
    d.insert(std::make_pair(peer, c));
  }

  typedef std::vector<int> VecInt;
  void print(const char* key, VecInt& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(VecInt, d, sItr)
      ss << *sItr << " ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  inline void setScalar(apf::Field* f, apf::MeshEntity* e, double val) {
    apf::setScalar(f,e,0,val);
  }
  inline double getScalar(apf::Field* f, apf::MeshEntity* e) {
    return apf::getScalar(f,e,0);
  }
  apf::Field* initScalarField(apf::Mesh* m, const char* name, 
      double initVal=0) {
    apf::Field* f = createFieldOn(m, name, apf::SCALAR);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while( (e = m->iterate(it)) )
      setScalar(f,e,initVal);
    m->end(it);
    return f;
  }
  apf::Field* initElmScalarField(apf::Mesh* m, const char* name, 
      double initVal=0) {
    const int dim = m->getDimension();
    apf::FieldShape* s = apf::getConstant(dim);
    apf::Field* f = createField(m, name, apf::SCALAR, s);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);
    while( (e = m->iterate(it)) )
      setScalar(f,e,initVal);
    m->end(it);
    return f;
  }

  inline void setNumber(apf::Numbering* n, apf::MeshEntity* e, 
      int val, int comp) {
    int node = 0;
    apf::number(n,e,node,comp,val);
  }
  inline int getNumber(apf::Numbering* n, apf::MeshEntity* e, int comp) {
    int node = 0;
    return apf::getNumber(n,e,node,comp);
  }
  apf::Numbering* initNumbering(apf::Mesh* m, const char* name, 
      int components, int initVal=0) {
    apf::FieldShape* s = m->getShape();
    apf::Numbering* n = apf::createNumbering(m,name,s,components);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while( (e = m->iterate(it)) )
      for(int c=0; c<components; c++)
        setNumber(n,e,initVal,c);
    m->end(it);
    return n;
  }

  apf::Numbering* initElmNumbering(apf::Mesh* m, const char* name, 
      int components, int initVal=0) {
    const int dim = m->getDimension();
    apf::FieldShape* s = apf::getConstant(dim);
    apf::Numbering* n = apf::createNumbering(m,name,s,components);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);
    while( (e = m->iterate(it)) )
      for(int c=0; c<components; c++)
        setNumber(n,e,initVal,c);
    m->end(it);
    return n;
  }
  void writeVtk(apf::Mesh* m, const char* pre) {
    std::stringstream ss;
    ss << pre << WeldSelectorCalls << "_";
    std::string s = ss.str();
    const int ids[] = {945,999,994,1051,827,829,831,811,814,951};
    const size_t len = sizeof(ids)/sizeof(int);
    int self = PCU_Comm_Self();
    for(size_t i=0; i<len; i++)
      if( ids[i] == self )
        apf::writeOneVtkFile(s.c_str(), m);
  }

  int const sdm[4] = {0,-1,1,1};
}

namespace parma {
  class WeldSelector : public Selector {
    public:
      WeldSelector(apf::Mesh* m, apf::MeshTag* wtag, Sides* sides)
        : Selector(m, wtag), s(sides)
      {
        WeldSelectorCalls++;
        int d = mesh->getDimension();
        int initVal = -1;
        nSrcDest[0] = initNumbering(mesh, "parmaVtx", 2, initVal);
        nSrcDest[sdm[d]] = 
          initElmNumbering(mesh, "parmaElm", 2, initVal);
        fSrcW[0] = initScalarField(mesh, "parmaSourceWeightVtx");
        fSrcW[sdm[d]] = initElmScalarField(mesh, "parmaSourceWeightElm");
      }
      ~WeldSelector() {
        apf::destroyNumbering(nSrcDest[0]);
        apf::destroyNumbering(nSrcDest[sdm[mesh->getDimension()]]);
        apf::destroyField(fSrcW[0]);
        apf::destroyField(fSrcW[sdm[mesh->getDimension()]]);
      }
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        select(tgts, plan);
        return plan;
      }
    protected:
      int getWeight(apf::MeshEntity* e) {
        int d = apf::Mesh::typeDimension[mesh->getType(e)];
        assert( d == mesh->getDimension() || d == 0 );
        return getScalar(fSrcW[sdm[d]],e);
      }
      void setWeight(apf::MeshEntity* e, double w) {
        int d = apf::Mesh::typeDimension[mesh->getType(e)];
        assert( d == mesh->getDimension() || d == 0 );
        setScalar(fSrcW[sdm[d]],e,w);
      }
      int getSrc(apf::MeshEntity* e) {
        int d = apf::Mesh::typeDimension[mesh->getType(e)];
        assert( d == mesh->getDimension() || d == 0 );
        return getNumber(nSrcDest[sdm[d]],e,0);
      }
      void setSrc(apf::MeshEntity* e, int src) {
        int d = apf::Mesh::typeDimension[mesh->getType(e)];
        assert( d == mesh->getDimension() || d == 0 );
        setNumber(nSrcDest[sdm[d]],e,src,0);
      }
      void setDest(apf::MeshEntity* e, int dest) {
        int d = apf::Mesh::typeDimension[mesh->getType(e)];
        assert( d == mesh->getDimension() || d == 0 );
        setNumber(nSrcDest[sdm[d]],e,dest,1);
      }
      int getDest(apf::MeshEntity* e) {
        int d = apf::Mesh::typeDimension[mesh->getType(e)];
        assert( d == mesh->getDimension() || d == 0 );
        return getNumber(nSrcDest[sdm[d]],e,0);
      }
      SetInt* getSmallSides(Targets* tgts, double selfW) {
        double maxW, avgW; 
        maxW = avgW = selfW;
        PCU_Max_Doubles(&maxW, 1);
        PCU_Add_Doubles(&avgW, 1);
        avgW /= PCU_Comm_Peers();
        const double tgtImb = 1 + (maxW/avgW - 1)*.75;
        const double tol = avgW * tgtImb;
        if( !PCU_Comm_Self() )
          fprintf(stdout, "DEBUG smallSide tolerance %.3f\n", tgtImb);

        int minSide = INT_MAX;
        if( selfW > tol ) {
          tgts->begin();
          const Targets::Item* t;
          while( (t = tgts->iterate()) ) {
            int peerSides = s->get(t->first);
            if( peerSides < minSide )
              minSide = peerSides;
          }
          tgts->end();
        }
        int small = minSide;
        PCU_Min_Ints(&minSide, 1);
        int minSideCount = ( small == minSide );
        PCU_Debug_Print("step %d small %d minSide %d tgts->size() %lu minSideCount %d\n",
            WeldSelectorCalls, small, minSide, tgts->size(), minSideCount);
        PCU_Add_Ints(&minSideCount, 1);
        if( !PCU_Comm_Self() )
          fprintf(stdout, "DEBUG minSide is %d minSide part count is %d\n", 
              minSide, minSideCount);

        SetInt* smallSidePeers = new SetInt;
        if( selfW > tol ) {
          s->begin();
          const Sides::Item* side;
          while( (side = s->iterate()) )
            if( side->second <= minSide )
              smallSidePeers->insert(side->first);
          s->end();
        }
        return smallSidePeers;
      }
      int getTarget(Targets* tgts, SetInt* smallSidePeers) {
        int maxPeer = -1;
        double maxW = -1;
        tgts->begin();
        const Targets::Item* t;
        while( (t = tgts->iterate()) ) {
          int peer = t->first;
          double peerW = t->second;
          if( smallSidePeers->count(peer) && peerW > maxW ) {
            maxPeer = peer;
            maxW = peerW;
          }
        }
        tgts->end();
        return maxPeer;
      }
      SetEnt* getJoint(int peer, double selfW) {
        SetEnt* joint = new SetEnt;
        apf::MeshEntity* v;
        apf::MeshIterator* it = mesh->begin(0);
        while( (v = mesh->iterate(it)) ) {
          apf::Copies rmts;
          mesh->getRemotes(v,rmts);
          if( rmts.count(peer) ) {
            joint->insert(v);
            setSrc(v,PCU_Comm_Self());
            setDest(v,peer);
            setWeight(v, selfW);
          }
        }
        mesh->end(it);
        return joint;
      }
      void select(Targets* tgts, apf::Migration* plan) {
        double selfW = parma::getWeight(mesh, wtag, mesh->getDimension());
        SetInt* smallSidePeers = getSmallSides(tgts,selfW);
        print("smallSidePeers", *smallSidePeers);
        int peer = getTarget(tgts,smallSidePeers);
        delete smallSidePeers;

        SetEnt* joint = getJoint(peer, selfW);
        writeVtk(mesh, "weldA");
        PCU_Comm_Begin();
        APF_ITERATE(SetEnt, *joint, u) {
          apf::Copies rmts;
          mesh->getRemotes(*u,rmts);
          int dest = getDest(*u);
          APF_ITERATE(apf::Copies, rmts, v) {
            int peer = v->first;
            apf::MeshEntity* rmt = v->second;
            PCU_COMM_PACK(peer, dest);
            PCU_COMM_PACK(peer, rmt);
            PCU_COMM_PACK(peer, selfW);
            PCU_Debug_Print("sending ent %p to %d srcW %.3f dest %d\n", 
                rmt, peer, selfW, dest);
          }
        }
        PCU_Comm_Send();
        while( PCU_Comm_Receive() ) {
          int src = PCU_Comm_Sender();
          int dest;
          apf::MeshEntity* v; 
          double srcW;
          PCU_COMM_UNPACK(dest);
          PCU_COMM_UNPACK(v);
          PCU_COMM_UNPACK(srcW);
          PCU_Debug_Print("received ent %p from %d srcW %.3f dest %d\n", 
              v, src, srcW, dest);
          joint->insert(v);
          if( srcW > getWeight(v) ) {
            setSrc(v, src);
            setDest(v, dest);
            setWeight(v, srcW);
          }
        }
        writeVtk(mesh, "weldB");

        SetEnt jointElms;
        const int self = PCU_Comm_Self();
        APF_ITERATE(SetEnt, *joint, u) {
          double srcW = getWeight(*u);
          int src = getSrc(*u);
          int dest = getDest(*u);
          if( self == src || self == dest ) continue;
          apf::Adjacent adjElms;
          mesh->getAdjacent(*u, mesh->getDimension(), adjElms);
          APF_ITERATE(apf::Adjacent, adjElms, elm) {
            if( srcW > getWeight(*elm) ) {
              setSrc(*elm, src);
              setDest(*elm, dest);
              setWeight(*elm, srcW);
              jointElms.insert(*elm);
            }
          }
        }
        writeVtk(mesh, "weldC");

        APF_ITERATE(SetEnt, jointElms, elm)
          plan->send(*elm, getDest(*elm));
        delete joint;
      }
    private:
      apf::Numbering* nSrcDest[2];
      apf::Field* fSrcW[2];
      WeldSelector();
      Weights* w;
      Sides* s;
  };

  Selector* makeWeldSelector(apf::Mesh* m, apf::MeshTag* wtag, Sides* s) {
    return new WeldSelector(m, wtag, s);
  }
}
