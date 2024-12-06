#include "parma_sides.h"
#include "parma_convert.h"
double parma::avgSharedSides(parma::Sides* s, pcu::PCU *PCUObj) {
  double tot = TO_DOUBLE(s->total());
  tot = PCUObj->Add<double>(tot);
  return tot / PCUObj->Peers();
}
