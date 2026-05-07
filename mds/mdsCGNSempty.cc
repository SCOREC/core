#include "apfMDS.h"

namespace apf {

Mesh2 *loadMdsFromCGNS(
  PCU_t, gmi_model*, const char*, apf::CGNSBCMap&,
  const std::vector<std::pair<std::string, std::string>>&
) {
  fail("Build with PUMI_ENABLE_CGNS to enable loadMdsFromCGNS.");
}

Mesh2 *loadMdsFromCGNS(
  PCU_t, gmi_model*, const char*, apf::CGNSBCMap&
) {
  fail("Build with PUMI_ENABLE_CGNS to enable loadMdsFromCGNS.");
}

} // namespace apf
