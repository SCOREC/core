#ifndef SAM_H
#define SAM_H

#include <apf.h>

namespace sam {

apf::Field* specifiedIso(apf::Mesh* m, const char* fieldName, const unsigned idx);
apf::Field* errorThreshold(apf::Mesh* m, const char* fieldName,
    const unsigned idx, const double limit, const double factor);

}

#endif
