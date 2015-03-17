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
    unsigned depth;
    double len;
  };
  bool compareComp(Comp a, Comp b) {
    if(a.depth > b.depth)
      return true;
    else 
      return false;
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
    for(unsigned i=0; i<n; i++) 
      depth[i] = 0;
    idT = m->createIntTag("parmaVtxCompId",1);
    markVertices();
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
    delete [] oldToNew;
  }

  void DCC::sortByDepth() {
    Comp* comp = new Comp[n];
    for(unsigned i=0; i<n; i++) {
      comp[i].i = i;
      comp[i].depth = depth[i];
    }
    std::stable_sort(comp, comp+n, compareComp);
    unsigned* order = new unsigned[n];
    for(unsigned i=0; i<n; i++)
      order[i] = comp[i].i;
    reorder(order);
    delete [] order;
    delete [] comp;
  }

  unsigned DCC::size() { return n; }

  unsigned DCC::iso() { return getNumIso(); }

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

  /**
   * brief get the core vertices in the component
   * remark if a component has no vertices assigned to it
   then it must has only boundary vertices
   set the core vertices to the bdry vertices in that case
   */
  void DCC::getCoreVerts() {
    for(unsigned i=0; i<size(); i++) {
      walkInward(i);
      if( !core[i].size() ) {
        PCU_Debug_Print("core %u is empty... assigning core to bdry\n", i);
        core[i] = bdry[i];
      }
    }
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
    for(unsigned i=0; i<size(); i++) {
      assert( core[i].size() );
      apf::MeshEntity* e = *(core[i].begin());
      core[i].clear();
      core[i].insert(e);
    }
  }

  void DCC::markVertices() {
    for(unsigned i=0; i<n; i++) {
      apf::MeshEntity* seed = getSeedEnt(i);
      walkComp(seed,i);
    }
  }

  /**
   * brief assign vertices to the component 
   * param src (In) element in the component
   * param comp (In) component id
   */
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
      apf::Downward verts;
      const int nv = m->getDownward(elm, 0, verts);
      setElmVtxIds(verts, nv, comp);
      addElmVtxToBdry(verts, nv, comp);
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
   * param verts (In) vertices to assign
   * param compId (In) id to assign to the vertices
   */
  void DCC::setElmVtxIds(apf::Downward& verts, const int nv, unsigned compId) {
    for(int i=0; i<nv; i++)
      if( !has(verts[i]) || (compId < getId(verts[i])) )
        setId(verts[i], compId);
  }

  void DCC::addElmVtxToBdry(apf::Downward& verts, const int nv, unsigned compId) {
    for(int i=0; i<nv; i++)
      if( onBoundary(m,verts[i]) )
        bdry[compId].insert(verts[i]);
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
  unsigned dcComponents::numIso() { return c->iso(); }
  bool dcComponents::has(apf::MeshEntity* e) { return c->has(e); }
  apf::MeshEntity* dcComponents::getCore(unsigned i) { return c->getCoreVtx(i); }
  bool dcComponents::bdryHas(unsigned i, apf::MeshEntity* e) {
    return c->getBdry(i)->count(e);
  }
  void dcComponents::beginBdry(unsigned i) { bItr->begin(c->getBdry(i)); }
  apf::MeshEntity* dcComponents::iterateBdry() { return bItr->iterate(); }
  void dcComponents::endBdry() { bItr->end(); }
} // end namespace parma
