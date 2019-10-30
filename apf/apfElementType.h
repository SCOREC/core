#ifndef APFELEMENTTYPE_H
#define APFELEMENTTYPE_H

#include "apfComplexType.h"

namespace apf
{

template <class T>
class ElementBase;
typedef ElementBase<double> Element;
typedef ElementBase<double_complex> ComplexElement;

};

#endif
