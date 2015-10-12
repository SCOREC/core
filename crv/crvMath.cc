/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crv.h"
#include "crvMath.h"
#include <mthQR.h>

namespace crv {

/*
 * compute inverse of A using QR to solve
 * Ai = Ri*Qt where Ri is obtained by back sub
 */

void invertMatrix(int n, mth::Matrix<double>& A,
    mth::Matrix<double>& Ai)
{
  mth::Matrix<double> Q(n,n), R(n,n), Ri(n,n);
  mth::decomposeQR(A,Q,R);

  mth::Vector<double> b(n);
  mth::Vector<double> x(n);
  for (int i = 0; i < n; ++i){
    b.zero();
    b[i] = 1.;
    mth::backsubUT(R,b,x);
    for (int j = 0; j < n; ++j)
      Ri(j,i) = x[j];
  }
  Ai.zero();

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        Ai(i,j) += Ri(i,k)*Q(j,k);
}

/*
 * 33 is the maximum in the table, before integer overflows occur
 * however, for n > 22, trinomial can overflow, and
 * For n > 18, quadnomial(n,i,j,k) can exceed MAX_INT and long's would be needed
 * This is also an upper bound on the order of tets, and implies a max order
 * of 7 to guarantee the full bezier jacobian determinant can work
 */
int binomial(int n, int i)
{

  i = std::min(n-i,i);

  if(i == 0)
    return 1;
  if(i == 1)
    return n;

  static int const bn4[1] = {6};
  static int const bn5[1] = {10};
  static int const bn6[2] = {15,20};
  static int const bn7[2] = {21,35};
  static int const bn8[3] = {28,56,70};
  static int const bn9[3] = {36,84,126};
  static int const bn10[4] = {45,120,210,252};
  static int const bn11[4] = {55,165,330,462};
  static int const bn12[5] = {66,220,495,792,924};
  static int const bn13[5] = {78,286,715,1287,1716};
  static int const bn14[6] = {91,364,1001,2002,3003,3432};
  static int const bn15[6] = {105,455,1365,3003,5005,6435};
  static int const bn16[7] = {120,560,1820,4368,8008,11440,12870};
  static int const bn17[7] = {136,680,2380,6188,12376,19448,24310};
  static int const bn18[8] = {153,816,3060,8568,18564,31824,43758,48620};
  static int const bn19[8] = {171,969,3876,11628,27132,50388,75582,92378};
  static int const bn20[9] = {190,1140,4845,15504,38760,77520,125970,167960,
      184756};
  static int const bn21[9] = {210,1330,5985,20349,54264,116280,203490,293930,
      352716};
  static int const bn22[10] = {231,1540,7315,26334,74613,170544,319770,497420,
      646646,705432};
  static int const bn23[10] = {253,1771,8855,33649,100947,245157,490314,817190,
      1144066,1352078};
  static int const bn24[11] = {276,2024,10626,42504,134596,346104,735471,1307504,
      1961256,2496144,2704156};
  static int const bn25[11] = {300,2300,12650,53130,177100,480700,1081575,2042975,
      3268760,4457400,5200300};
  static int const bn26[12] = {325,2600,14950,65780,230230,657800,1562275,3124550,
      5311735,7726160,9657700,10400600};
  static int const bn27[12] = {351,2925,17550,80730,296010,888030,2220075,4686825,
      8436285,13037895,17383860,20058300};
  static int const bn28[13] = {378,3276,20475,98280,376740,1184040,3108105,6906900,
      13123110,21474180,30421755,37442160,40116600};
  static int const bn29[13] = {406,3654,23751,118755,475020,1560780,4292145,10015005,
      20030010,34597290,51895935,67863915,77558760};
  static int const bn30[14] = {435,4060,27405,142506,593775,2035800,5852925,14307150,
      30045015,54627300,86493225,119759850,145422675,155117520};
  static int const bn31[14] = {465,4495,31465,169911,736281,2629575,7888725,20160075,
      44352165,84672315,141120525,206253075,265182525,300540195};
  static int const bn32[15] = {496,4960,35960,201376,906192,3365856,10518300,28048800,
      64512240,129024480,225792840,347373600,471435600,565722720,601080390};
  static int const bn33[15] = {528,5456,40920,237336,1107568,4272048,13884156,38567100,
      92561040,193536720,354817320,573166440,818809200,1037158320,1166803110};


  static int const* const bnTable[34] = {0,0,0,0,bn4,bn5,bn6,bn7,bn8,
      bn9,bn10,bn11,bn12,bn13,bn14,bn15,bn16,bn17,bn18,bn19,bn20,bn21,bn22,bn23,
      bn24,bn25,bn26,bn27,bn28,bn29,bn30,bn31,bn32,bn33};

  return bnTable[n][i-2];

}

int trinomial(int n, int i, int j)
{
  return binomial(n,i)*binomial(n-i,j);
}

int quadnomial(int n, int i, int j, int k)
{
  return binomial(n,i)*binomial(n-i,j)*binomial(n-i-j,k);
}

}
