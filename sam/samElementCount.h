#ifndef SAM_ELEMENT_COUNT_H
#define SAM_ELEMENT_COUNT_H

namespace apf {
  class Field;
}

namespace sam {

double getIsoLengthScalar(apf::Field* iso_field, double targetElementCount);
void scaleIsoSizeField(apf::Field* iso_field, double targetElementCount);

}

#endif
