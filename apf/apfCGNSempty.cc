#include <apf.h>

namespace apf {

void writeCGNS(const char*, Mesh*, const apf::CGNSBCMap&) {
  fail("Build with PUMI_ENABLE_CGNS to enable apf::writeCGNS.");
}

} // namespace apf
