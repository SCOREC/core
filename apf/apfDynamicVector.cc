#include "apfDynamicVector.h"

std::ostream& operator<<(std::ostream& s, apf::DynamicVector const& x)
{
  for (std::size_t i = 0; i < x.getSize(); ++i) {
      s << x(i) << '\n';
  }
  return s;
}
