#include <parma.h>
#include <parma_ghostOwner.h>
#include <sstream>
#include <fstream>
#include <pcu_util.h>

namespace {
  class Ptn {
    private:
    int c;
    int pp;
    int f;
    int countOwned(apf::Mesh* m) {
      int owned = 0;
      apf::MeshEntity* vtx;
      apf::MeshIterator* itr = m->begin(0);
      while( (vtx = m->iterate(itr)) )
        owned += parma::isOwned(m,vtx);
      m->end(itr);
      return owned;
    }
    public:
    Ptn(apf::Mesh* m) {
      const long totv = m->getPCU()->Add<long>(countOwned(m));
      c = pp = totv / m->getPCU()->Peers();
      f = pp * m->getPCU()->Self();
      const int remainder = totv % m->getPCU()->Peers();
      if( m->getPCU()->Self() == m->getPCU()->Peers()-1 )
        c += remainder;
    }
    int getWriter(int id, pcu::PCU *PCUObj) {
      int writer = id / pp;
      if ( writer == PCUObj->Peers() )
        writer--;
      return writer;
    }
    int count() {
      return c;
    }
    int first() {
      return f;
    }
  };

  void getPtnArray(apf::Mesh* m, Ptn& p, int* ptn) {
    apf::MeshTag* t = m->findTag("ugrid-vtx-ids");
    PCU_ALWAYS_ASSERT(t);
    apf::MeshEntity* vtx;
    apf::MeshIterator* itr = m->begin(0);
    int id = 0;
    m->getPCU()->Begin();
    while( (vtx = m->iterate(itr)) ) {
      if( parma::isOwned(m, vtx) ) {
        m->getIntTag(vtx, t, &id);
        m->getPCU()->Pack(p.getWriter(id, m->getPCU()), id);
      }
    }
    m->end(itr);
    m->getPCU()->Send();
    while( m->getPCU()->Receive() ) {
      int id = 0;
      m->getPCU()->Unpack(id);
      const int idx = id - p.first();
      PCU_ALWAYS_ASSERT(idx >= 0 && idx < p.count());
      ptn[idx] = m->getPCU()->Sender();
    }
  }

  void writePtnArray(int* a, int n, std::fstream& f) {
    for(int i=0; i<n; i++)
      f << a[i] << '\n';
  }

  void open(const char* name, std::fstream& f, pcu::PCU *PCUObj) {
    std::stringstream ss;
    ss << name << PCUObj->Self() << ".ptn";
    std::string s = ss.str();
    f.open(s.c_str(), std::fstream::out);
  }

  void writeVtxPtn(apf::Mesh* m, const char* name) {
    PCU_ALWAYS_ASSERT(name);
    std::fstream f;
    open(name,f,m->getPCU());
    Ptn p(m);
    int* ptn = new int[p.count()];
    getPtnArray(m, p, ptn);
    writePtnArray(ptn, p.count(), f);
    f.close();
    delete [] ptn;
  }
} //end namespace

void Parma_WriteVtxPtn(apf::Mesh* m, const char* prefix) {
  writeVtxPtn(m,prefix);
}
