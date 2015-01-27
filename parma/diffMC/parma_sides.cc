#include <PCU.h>
#include "parma_sides.h"
double parma::avgSharedSides(parma::Sides* s) {
  double tot = static_cast<double>(s->total());
  PCU_Add_Doubles(&tot, 1);
  return tot / PCU_Comm_Peers();
}
