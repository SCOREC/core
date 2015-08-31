/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvTables.h"

namespace crv {

/* beyond this comment lie the vast tables
 * in which the Programmer in His wisdom
 * learned from the master Programmer,
 * who in His wisdom, laid down this construct for tables,
 * repurposed to encode lookup tables for bezier tris,
 * which are formulated as [P][i+j*(P+1)-j*(j-1)/2],
 * with k implied, because I was not yet learned at the time,
 * and a lookup table for bezier tets,
 * as a [P][i][j][k] formulation, with l implied,
 * and while first order is superfluous, consistency is maintained
 */

// Notes: 0 is a fake point, added for de Casteljau's algorithm
static unsigned const b2_0_0[1] = {2};
static unsigned const* const b2_0[1] = {b2_0_0};
static unsigned const b2_1_0[2] = {2,1};
static unsigned const b2_1_1[1] = {0};
static unsigned const* const b2_1[2] = {b2_1_0,b2_1_1};
static unsigned const b2_2_0[3] = {2,4,1};
static unsigned const b2_2_1[2] = {5,3};
static unsigned const b2_2_2[1] = {0};
static unsigned const* const b2_2[3] = {b2_2_0,b2_2_1,b2_2_2};
static unsigned const b2_3_0[4] = {2,6,5,1};
static unsigned const b2_3_1[3] = {7,9,4};
static unsigned const b2_3_2[2] = {8,3};
static unsigned const b2_3_3[1] = {0};
static unsigned const* const b2_3[4] = {b2_3_0,b2_3_1,b2_3_2,b2_3_3};
static unsigned const b2_4_0[5] = {2,8,7,6,1};
static unsigned const b2_4_1[4] = {9,14,13,5};
static unsigned const b2_4_2[3] = {10,12,4};
static unsigned const b2_4_3[2] = {11,3};
static unsigned const b2_4_4[1] = {0};
static unsigned const* const b2_4[5] = {b2_4_0,b2_4_1,b2_4_2,b2_4_3,b2_4_4};
static unsigned const b2_5_0[6] = {2,10,9,8,7,1};
static unsigned const b2_5_1[5] = {11,19,18,17,6};
static unsigned const b2_5_2[4] = {12,20,16,5};
static unsigned const b2_5_3[3] = {13,15,4};
static unsigned const b2_5_4[2] = {14,3};
static unsigned const b2_5_5[1] = {0};
static unsigned const* const b2_5[6] =
{b2_5_0,b2_5_1,b2_5_2,b2_5_3,b2_5_4,b2_5_5};
static unsigned const b2_6_0[7] = {2,12,11,10,9,8,1};
static unsigned const b2_6_1[6] = {13,24,23,22,21,7};
static unsigned const b2_6_2[5] = {14,25,27,20,6};
static unsigned const b2_6_3[4] = {15,26,19,5};
static unsigned const b2_6_4[3] = {16,18,4};
static unsigned const b2_6_5[2] = {17,3};
static unsigned const b2_6_6[1] = {0};
static unsigned const* const b2_6[7] =
{b2_6_0,b2_6_1,b2_6_2,b2_6_3,b2_6_4,b2_6_5,b2_6_6};
static unsigned const b2_7_0[8] = {2,14,13,12,11,10,9,1};
static unsigned const b2_7_1[7] = {15,29,28,27,26,25,8};
static unsigned const b2_7_2[6] = {16,30,35,34,24,7};
static unsigned const b2_7_3[5] = {17,31,33,23,6};
static unsigned const b2_7_4[4] = {18,32,22,5};
static unsigned const b2_7_5[3] = {19,21,4};
static unsigned const b2_7_6[2] = {20,3};
static unsigned const b2_7_7[1] = {0};
static unsigned const* const b2_7[8] =
{b2_7_0,b2_7_1,b2_7_2,b2_7_3,b2_7_4,b2_7_5,b2_7_6,b2_7_7};
static unsigned const b2_8_0[9] = {2,16,15,14,13,12,11,10,1};
static unsigned const b2_8_1[8] = {17,34,33,32,31,30,29,9};
static unsigned const b2_8_2[7] = {18,35,43,42,41,28,8};
static unsigned const b2_8_3[6] = {19,36,44,40,27,7};
static unsigned const b2_8_4[5] = {20,37,39,26,6};
static unsigned const b2_8_5[4] = {21,38,25,5};
static unsigned const b2_8_6[3] = {22,24,4};
static unsigned const b2_8_7[2] = {23,3};
static unsigned const b2_8_8[1] = {0};
static unsigned const* const b2_8[9] =
{b2_8_0,b2_8_1,b2_8_2,b2_8_3,b2_8_4,b2_8_5,b2_8_6,b2_8_7,b2_8_8};
static unsigned const b2_9_0[10] = {2,18,17,16,15,14,13,12,11,1};
static unsigned const b2_9_1[9] = {19,39,38,37,36,35,34,33,10};
static unsigned const b2_9_2[8] = {20,40,51,50,49,48,32,9};
static unsigned const b2_9_3[7] = {21,41,52,54,47,31,8};
static unsigned const b2_9_4[6] = {22,42,53,46,30,7};
static unsigned const b2_9_5[5] = {23,43,45,29,6};
static unsigned const b2_9_6[4] = {24,44,28,5};
static unsigned const b2_9_7[3] = {25,27,4};
static unsigned const b2_9_8[2] = {26,3};
static unsigned const b2_9_9[1] = {0};
static unsigned const* const b2_9[10] =
{b2_9_0,b2_9_1,b2_9_2,b2_9_3,b2_9_4,b2_9_5,b2_9_6,b2_9_7,b2_9_8,b2_9_9};
static unsigned const b2_10_0[11] = {2,20,19,18,17,16,15,14,13,12,1};
static unsigned const b2_10_1[10] = {21,44,43,42,41,40,39,38,37,11};
static unsigned const b2_10_2[9] = {22,45,59,58,57,56,55,36,10};
static unsigned const b2_10_3[8] = {23,46,60,65,64,54,35,9};
static unsigned const b2_10_4[7] = {24,47,61,63,53,34,8};
static unsigned const b2_10_5[6] = {25,48,62,52,33,7};
static unsigned const b2_10_6[5] = {26,49,51,32,6};
static unsigned const b2_10_7[4] = {27,50,31,5};
static unsigned const b2_10_8[3] = {28,30,4};
static unsigned const b2_10_9[2] = {29,3};
static unsigned const b2_10_10[1] = {0};
static unsigned const* const b2_10[11] =
{b2_10_0,b2_10_1,b2_10_2,b2_10_3,b2_10_4,b2_10_5,b2_10_6,b2_10_7,b2_10_8,
    b2_10_9,b2_10_10};

unsigned const* const* const b2[11] =
{b2_0,b2_1,b2_2,b2_3,b2_4,b2_5,b2_6,b2_7,b2_8,b2_9,b2_10};

static unsigned const b3_0_00[1] = {3};
static unsigned const* const b3_0_0[1] = {b3_0_00};
static unsigned const* const* const b3_0[1] = {b3_0_0};

static unsigned const b3_1_00[2] = {3,2};
static unsigned const b3_1_01[1] = {1};
static unsigned const b3_1_10[1] = {0};
static unsigned const* const b3_1_0[2] = {b3_1_00,b3_1_01};
static unsigned const* const b3_1_1[1] = {b3_1_10};
static unsigned const* const* const b3_1[2] = {b3_1_0,b3_1_1};

static unsigned const b3_2_00[3] = {3,9,2};
static unsigned const b3_2_01[2] = {8,5};
static unsigned const b3_2_02[1] = {1};
static unsigned const b3_2_10[2] = {7,6};
static unsigned const b3_2_11[1] = {4};
static unsigned const b3_2_20[1] = {0};
static unsigned const* const b3_2_0[3] = {b3_2_00,b3_2_01,b3_2_02};
static unsigned const* const b3_2_1[2] = {b3_2_10,b3_2_11};
static unsigned const* const b3_2_2[1] = {b3_2_20};
static unsigned const* const* const b3_2[3] = {b3_2_0,b3_2_1,b3_2_2};

static unsigned const b3_3_00[4] = {3,15,14,2};
static unsigned const b3_3_01[3] = {13,18,7};
static unsigned const b3_3_02[2] = {12,6};
static unsigned const b3_3_03[1] = {1};
static unsigned const b3_3_10[3] = {11,19,8};
static unsigned const b3_3_11[2] = {17,16};
static unsigned const b3_3_12[1] = {5};
static unsigned const b3_3_20[2] = {10,9};
static unsigned const b3_3_21[1] = {4};
static unsigned const b3_3_30[1] = {0};
static unsigned const* const b3_3_0[4] = {b3_3_00,b3_3_01,b3_3_02,b3_3_03};
static unsigned const* const b3_3_1[3] = {b3_3_10,b3_3_11,b3_3_12};
static unsigned const* const b3_3_2[2] = {b3_3_20,b3_3_21};
static unsigned const* const b3_3_3[1] = {b3_3_30};
static unsigned const* const* const b3_3[4] =
  {b3_3_0,b3_3_1,b3_3_2,b3_3_3};

static unsigned const b3_4_00[5] = {3,21,20,19,2};
static unsigned const b3_4_01[4] = {18,30,29,9};
static unsigned const b3_4_02[3] = {17,28,8};
static unsigned const b3_4_03[2] = {16,7};
static unsigned const b3_4_04[1] = {1};
static unsigned const b3_4_10[4] = {15,33,32,10};
static unsigned const b3_4_11[3] = {27,34,24};
static unsigned const b3_4_12[2] = {26,23};
static unsigned const b3_4_13[1] = {6};
static unsigned const b3_4_20[3] = {14,31,11};
static unsigned const b3_4_21[2] = {25,22};
static unsigned const b3_4_22[1] = {5};
static unsigned const b3_4_30[2] = {13,12};
static unsigned const b3_4_31[1] = {4};
static unsigned const b3_4_40[1] = {0};
static unsigned const* const b3_4_0[5] =
  {b3_4_00,b3_4_01,b3_4_02,b3_4_03,b3_4_04};
static unsigned const* const b3_4_1[4] = {b3_4_10,b3_4_11,b3_4_12,b3_4_13};
static unsigned const* const b3_4_2[3] = {b3_4_20,b3_4_21,b3_4_22};
static unsigned const* const b3_4_3[2] = {b3_4_30,b3_4_31};
static unsigned const* const b3_4_4[1] = {b3_4_40};
static unsigned const* const* const b3_4[5] =
  {b3_4_0,b3_4_1,b3_4_2,b3_4_3,b3_4_4};

unsigned const* const* const* const b3[5] =
{b3_0,b3_1,b3_2,b3_3,b3_4};

// publically accessible access
int getTriPointIndex(int P, int i, int j)
{
  return crv::b2[P][i][j];
}

int getTetPointIndex(int P, int i, int j, int k)
{
  return crv::b3[P][i][j][k];
}

}
