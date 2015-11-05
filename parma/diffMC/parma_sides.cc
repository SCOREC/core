#include <PCU.h>
#include "parma_sides.h"
double parma::avgSharedSides(parma::Sides* s) {
  double tot = static_cast<double>(s->total());
  tot = PCU_Add_Double(tot);
  return tot / PCU_Comm_Peers();
}
