#include <PCU.h>
#include "parma_sides.h"
#include "parma_convert.h"
double parma::avgSharedSides(parma::Sides* s) {
  double tot = TO_DOUBLE(s->total());
  tot = PCU_Add_Double(tot);
  return tot / PCU_Comm_Peers();
}
