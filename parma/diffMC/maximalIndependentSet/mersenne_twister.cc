#include "mersenne_twister.h"
#include "PCU.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Matsumoto, Makoto, and Takuji Nishimura.
   "Mersenne twister: a 623-dimensionally equidistributed
    uniform pseudo-random number generator."
   ACM Transactions on Modeling and Computer Simulation
   (TOMACS) 8.1 (1998): 3-30. */

/* borrowed from
   https://raw.githubusercontent.com/ibaned/tetknife/master/mersenne_twister.c
*/

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7FFFFFFF /* less significant w-r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)

static unsigned long mt[N]; /* the array for the state vector */
static int mti = N + 1; /* mti == N + 1 means mt is not initialized */

namespace {
  void fail(const char* msg) {
    if( ! PCU_Comm_Self() )
      fprintf(stderr, "%s", msg);
    exit(EXIT_FAILURE);
  }
}

void mersenne_twister_seed(unsigned seed)
{
  assert(seed);
  /* setting initial seeds to mt[N] using
     the generator Line 25 of Table 1 in
     [KNUTH 1981, The Art of Computer Programming
        Vol. 2 (2nd Ed), pp102] */
  mt[0] = seed;
  for (mti = 1; mti < N; ++mti)
    mt[mti] = (6909 * mt[mti - 1]) & 0xFFFFFFFF;
}

unsigned mersenne_twister(void)
{
  static unsigned long const mag01[2] = {0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A fox x = {0,1} */
  unsigned long y;
  int kk;
  if (mti == N + 1)
    fail("mersenne twister was not seeded before use\n");
  if (mti == N) { /* generate N words at one time */
    for (kk = 0; kk < N - M; ++kk) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 1];
    }
    for (; kk < N - 1; ++kk) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1];
    }
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 1];
    mti = 0;
  }
  y = mt[mti++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y ^= TEMPERING_SHIFT_L(y);
  return (unsigned) y;
}
