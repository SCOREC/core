#include "parma_selector.h"
#include "parma_targets.h"

namespace parma {
  class HpsSelector : public Selector {
    public:
      HpsSelector(apf::Mesh* m, apf::MeshTag* w)
        : Selector(m, w) {}
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        //compute the optimal threshold for merging and splitting
        //merge and split using computed threshold
        return plan;
      }
    private:
      apf::Mesh* mesh;
      apf::MeshTag* wtag;
      HpsSelector();
  };
  Selector* makeHpsSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new HpsSelector(m, w);
  }
} //end namespace parma
