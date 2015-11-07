#include "phKitchen.h"
#include "ph.h"
#include <ma.h>
#include <PCU.h>

namespace kitchen {
  void readAndAttachFields(gmi_model*&, apf::Mesh2*&, ph::Input&) {
  }
  void readAndAttachFields(gmi_model*&, apf::Mesh2*&,
      ph::Input&, RStream*) {
  }
  void preprocess(gmi_model*&, apf::Mesh2*&, ph::Input&, GRStream*) {
  }

  void adapt(gmi_model*, apf::Mesh2*& m, ma::IsotropicFunction* szFld) {
    ma::Input* ma_in = ma::configure(m, szFld);
    if (m->hasMatching()) {
      if (!PCU_Comm_Self())
        printf("Matched mesh: disabling coarsening, snapping, and shape correction,\n"
            "  synchronizing \"errors\" field (source of size)\n");
      ma_in->shouldCoarsen = false;
      ma_in->shouldSnap = false;
      ma_in->shouldFixShape = false;
    }
    ma::adapt(ma_in);
  }
}
