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
      const long totv = PCU_Add_Long(countOwned(m));
      c = pp = totv / PCU_Comm_Peers();
      f = pp * PCU_Comm_Self();
      const int remainder = totv % PCU_Comm_Peers();
      if( PCU_Comm_Self() == PCU_Comm_Peers()-1 )
        c += remainder;
    }
    int getWriter(int id) {
      int writer = id / pp;
      if ( writer == PCU_Comm_Peers() )
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
    PCU_Comm_Begin();
    while( (vtx = m->iterate(itr)) ) {
      if( parma::isOwned(m, vtx) ) {
        m->getIntTag(vtx, t, &id);
        PCU_COMM_PACK(p.getWriter(id), id);
      }
    }
    m->end(itr);
    PCU_Comm_Send();
    while( PCU_Comm_Receive() ) {
      int id = 0;
      PCU_COMM_UNPACK(id);
      const int idx = id - p.first();
      PCU_ALWAYS_ASSERT(idx >= 0 && idx < p.count());
      ptn[idx] = PCU_Comm_Sender();
    }
  }

  void writePtnArray(int* a, int n, std::fstream& f) {
    for(int i=0; i<n; i++)
      f << a[i] << '\n';
  }

  void open(const char* name, std::fstream& f) {
    std::stringstream ss;
    ss << name << PCU_Comm_Self() << ".ptn";
    std::string s = ss.str();
    f.open(s.c_str(), std::fstream::out);
  }

  void writeVtxPtn(apf::Mesh* m, const char* name) {
    PCU_ALWAYS_ASSERT(name);
    std::fstream f;
    open(name,f);
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
