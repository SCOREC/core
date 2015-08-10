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

static unsigned const b2_1[3] = {2,0,1};
static unsigned const b2_2[6] = {2,5,0,4,3,1};
static unsigned const b2_3[10] = {2,7,8,0,6,9,3,5,4,1};
static unsigned const b2_4[15] = {2,9,10,11,0,8,14,12,3,7,13,4,6,5,1};
static unsigned const b2_5[21] = {2,11,12,13,14,0,10,19,20,15,3,9,18,
  16,4,8,17,5,7,6,1};
static unsigned const b2_6[28] = {2,13,14,15,16,17,0,12,24,25,26,18,3,
  11,23,27,19,4,10,22,20,5,9,21,6,8,7,1};
unsigned const* const b2[7] =
{0,b2_1,b2_2,b2_3,b2_4,b2_5,b2_6};


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
{0,b3_1,b3_2,b3_3,b3_4};

}
