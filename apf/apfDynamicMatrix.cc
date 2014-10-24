#include "apfDynamicMatrix.h"

std::ostream& operator<<(std::ostream& s, apf::DynamicMatrix const& A)
{
  for (std::size_t i = 0; i < A.getRows(); ++i) {
    for (std::size_t j = 0; j < A.getColumns(); ++j)
      s << A(i,j) << ' ';
    s << '\n';
  }
  return s;
}
