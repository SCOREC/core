#ifndef SAM_H
#define SAM_H

#include <apf.h>

namespace sam {

apf::Field* specifiedIso(apf::Mesh* m, const char* fieldName, const unsigned idx);
apf::Field* compareIsoSF(apf::Mesh* m, const char* desiredSzFld, int method = 0);
apf::Field* errorThreshold(apf::Mesh* m, const char* fieldName,
    const unsigned idx, const double limit, const double factor);
void multiplySF(apf::Mesh* m, apf::Field* sf, double factor);
void multiplySFBox(apf::Mesh* m, apf::Field* sf, double factor, double* box);
void multiplySFCyl(apf::Mesh* m, apf::Field* sf, double factor, double* cyl);
void multiplySFRegion(apf::Mesh* m, apf::Field* sf, double factor, int tag);

}

#endif
