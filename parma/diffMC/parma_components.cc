#include <list>
#include <set>
#include <apf.h>
#include <PCU.h>
#include "parma_components.h"
#include "parma_meshaux.h"

#define TO_UINT(a) static_cast<unsigned>(a)
#define TO_INT(a) static_cast<int>(a)

namespace {
  struct Comp {
    unsigned i;
    unsigned d;
  };
  bool compareComp(Comp a, Comp b) {
    return (a.d > b.d);
  }

  void reduce(apf::Mesh* m, parma::Level* core) {
    double max = 0;
    apf::MeshEntity* maxVtx = NULL;
    APF_ITERATE(parma::Level, *core, e) {
      apf::Vector3 u;
      m->getPoint(*e, 0, u);
      double len = u.getLength();
      if( len > max ) {
        max = len;
        maxVtx = *e;
      }
    }
    assert(maxVtx);
    core->clear();
    core->insert(maxVtx);
  }
}

#define DCC dcComponents::Components
namespace parma {
  DCC::Components(apf::Mesh* mesh, unsigned verbose) 
    : dcPart(mesh,verbose), m(mesh), vb(verbose) 
  {
    n = getNumComps();
    depth = new unsigned[n];
    bdry = new Level[n];
    core = new Level[n];
    for(unsigned i=0; i<n; i++) depth[i] = 0;
    idT = m->createIntTag("parmaVtxCompId",1);
    markVertices();
    getBdryVerts();
    getCoreVerts();
    getCoreVtx();
    sortByDepth();
  }

  DCC::~Components() {
    delete [] depth;
    delete [] bdry;
    delete [] core;
    apf::removeTagFromDimension(m,idT,0);
    m->destroyTag(idT);
  }

  void DCC::reorder(unsigned* order) {
    unsigned* oldToNew = new unsigned[n];
    unsigned* dtmp = new unsigned[n];
    Level* btmp = new Level[n];
    apf::MeshEntity** ctmp = new apf::MeshEntity*[n];
    for(unsigned i=0; i<n; i++) {
      oldToNew[order[i]] = i;
      dtmp[i] = depth[i];
      btmp[i] = bdry[i];
      assert(1 == core[i].size());
      ctmp[i] = *(core[i].begin());
    }
    for(unsigned i=0; i<n; i++) {
      depth[i] = dtmp[order[i]];
      bdry[i] = btmp[order[i]];
      core[i].clear();
      core[i].insert(ctmp[order[i]]);
    }
    delete [] dtmp;
    delete [] btmp;
    delete [] ctmp;

    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) )
      if( has(v) ) {
        unsigned cid = getId(v);
        setId(v,oldToNew[cid]);
      }
    m->end(it);
  }

  void DCC::sortByDepth() {
    Comp* comp = new Comp[n];
    for(unsigned i=0; i<n; i++) {
      comp[i].i = i;
      comp[i].d = depth[i];
    }
    std::sort(comp, comp+n, compareComp);
    unsigned* order = new unsigned[n];
    for(unsigned i=0; i<n; i++) 
      order[i] = comp[i].i;
    reorder(order);
  }

  unsigned DCC::size() { return n; }

  unsigned DCC::getDepth(unsigned i) { assert(i<n); return depth[i]; }

  void DCC::setDepth(unsigned i, unsigned d) { assert(i<n); depth[i] = d; }

  Level* DCC::getBdry(unsigned i) { assert(i<n); return &(bdry[i]); }

  Level* DCC::getCore(unsigned i) { assert(i<n); return &(core[i]); }

  bool DCC::has(apf::MeshEntity* e) { return m->hasTag(e, idT); }

  apf::MeshEntity* DCC::getCoreVtx(unsigned i) { 
    assert(i<n); 
    Level* lvl = getCore(i);
    assert(1 == lvl->size());
    return *(lvl->begin());
  }

  unsigned DCC::getId(apf::MeshEntity* e) { 
    assert(m->hasTag(e, idT));
    int ctv; m->getIntTag(e, idT, &ctv);
    unsigned cid = TO_UINT(ctv);
    assert(cid < n);
    return cid;
  }

  void DCC::setId(apf::MeshEntity* e, unsigned compid) { 
    int cid = TO_INT(compid);
    m->setIntTag(e, idT, &cid);
  }

  void DCC::getBdryVerts() {
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while( (e = m->iterate(it)) )
      if( onBoundary(m,e) && has(e) ) {
        unsigned cid = getId(e);
        parma::Level* b = getBdry(cid);
        b->insert(e);
      }
    m->end(it);
  }

  void DCC::getCoreVerts() {
    for(unsigned i=0; i<size(); i++)
      walkInward(i);
  }

  /**
   * brief starting from the boundary vertices of the given component
   *       run a breadth first search
   * remark (0) set the 'core' level in the component to the leave nodes
   *        in the walk tree
   *        (1) set the 'depth' field of the component to tree's height
   * parma compId (in) component id of interest
   */
  void DCC::walkInward(unsigned compId) {
    apf::MeshTag* lvlT = m->createIntTag("parmaWalkLevels",1);
    parma::Level cur;
    parma::Level next = *(getBdry(compId));
    int treeDepth = 0;
    while( ! next.empty() ) {
      cur = next;
      next.clear();
      treeDepth++;
      APF_ITERATE(parma::Level, cur, vtxItr) {
        apf::MeshEntity* v = *vtxItr;
        if( m->hasTag(v,lvlT) ) continue;
        m->setIntTag(v,lvlT,&treeDepth);
        apf::Adjacent adjVtx;
        getEdgeAdjVtx(m,v,adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, vItr)
          if( ! m->hasTag(*vItr,lvlT) && 
              ! cur.count(*vItr) &&
              has(*vItr) &&
              (compId == getId(*vItr)) )
            next.insert(*vItr);
      }
    }
    parma::Level* lvl = getCore(compId);
    APF_ITERATE(parma::Level, cur, lItr) 
      lvl->insert(*lItr); // (0)
    setDepth(compId, TO_UINT(treeDepth)); // (1)

    apf::removeTagFromDimension(m,lvlT,0);
    m->destroyTag(lvlT);
  }

  void DCC::getCoreVtx() {
    for(unsigned i=0; i<size(); i++)
      reduce(m,getCore(i));
  }

  void DCC::markVertices() {
    for(unsigned i=0; i<n; i++) {
      apf::MeshEntity* seed = getSeedEnt(i);
      walkComp(seed,i);
    }
  }

  void DCC::walkComp(apf::MeshEntity* src, unsigned comp) {
    int one = 1;
    apf::MeshTag* vtag = m->createIntTag("walkCompVisited",1);

    std::list<apf::MeshEntity*> elms;
    elms.push_back(src);
    while( ! elms.empty() ) {
      apf::MeshEntity* elm = elms.front();
      elms.pop_front();
      if( m->hasTag(elm, vtag) ) continue;
      m->setIntTag(elm, vtag, &one);
      setElmVtxIds(elm, comp);
      apf::Adjacent adjElms;
      getDwn2ndAdj(m, elm, adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, eit)
        if( ! isIsolated(*eit) &&
            (compId(*eit) == comp) &&
            ! m->hasTag(*eit,vtag) )
          elms.push_back(*eit);
    }

    apf::removeTagFromDimension(m,vtag,m->getDimension());
    m->destroyTag(vtag);
  }

  /**
   * brief set vtx component ids
   * remark vertices shared by multiple components are assigned to the 
   *        component with the lowest id
   * param elm (In) elm with vertices to assign
   * param compId (In) id to assign to the vertices
   */
  void DCC::setElmVtxIds(apf::MeshEntity* elm, unsigned compId) {
    apf::Downward vtx;
    int nv = m->getDownward(elm, 0, vtx);
    for(int i=0; i<nv; i++)
      if( !has(vtx[i]) || (compId < getId(vtx[i])) )
        setId(vtx[i], compId);
  }

  class dcComponents::BdryItr {
    public:
      BdryItr() : active(false) {}
      void begin(Level* l) {
        assert(!active);
        lvl = l;
        active = true;
        itr = lvl->begin();
      }
      apf::MeshEntity* iterate() {
        assert(active);
        if( itr == lvl->end() ) 
          return NULL;
        else
          return *(itr++);
      }
      void end() {
        assert(active);
        active = false;
      }
    private:
      bool active;
      Level::iterator itr;
      Level* lvl;
  };

  dcComponents::dcComponents(apf::Mesh* m, unsigned verbose) 
    : c(new Components(m,verbose)), bItr(new BdryItr) {}
  dcComponents::~dcComponents() { delete c; delete bItr; }
  unsigned dcComponents::size() { return c->size(); }
  unsigned dcComponents::getId(apf::MeshEntity* e) { return c->getId(e); }
  bool dcComponents::has(apf::MeshEntity* e) { return c->has(e); }
  apf::MeshEntity* dcComponents::getCore(unsigned i) { return c->getCoreVtx(i); }
  void dcComponents::beginBdry(unsigned i) { bItr->begin(c->getBdry(i)); }
  apf::MeshEntity* dcComponents::iterateBdry() { return bItr->iterate(); }
  void dcComponents::endBdry() { bItr->end(); }
} // end namespace parma
