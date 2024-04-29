#ifndef MERSENNE_TWISTER_H
#define MERSENNE_TWISTER_H

#include <PCUObj.h>

void mersenne_twister_seed(unsigned seed);
unsigned mersenne_twister(pcu::PCU *PCUObj);

#endif
