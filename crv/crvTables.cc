/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvTables.h"
#include "apf.h"
#include "apfMesh.h"

namespace crv {

/* beyond this comment lie the vast tables
 * in which the Programmer in His wisdom
 * learned from the master Programmer,
 * who in His wisdom, laid down this construct for tables,
 * repurposed to encode lookup tables for bezier tris,
 * which are formulated as [P][i][j],
 * with k implied, because I was not yet learned at the time,
 * and a lookup table for bezier tets,
 * as a [P][i][j][k] formulation, with l implied,
 * and while zeroth and first order are superfluous, consistency is maintained
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
static unsigned const b2_5_1[5] = {11,20,19,17,6};
static unsigned const b2_5_2[4] = {12,18,16,5};
static unsigned const b2_5_3[3] = {13,15,4};
static unsigned const b2_5_4[2] = {14,3};
static unsigned const b2_5_5[1] = {0};
static unsigned const* const b2_5[6] =
{b2_5_0,b2_5_1,b2_5_2,b2_5_3,b2_5_4,b2_5_5};
static unsigned const b2_6_0[7] = {2,12,11,10,9,8,1};
static unsigned const b2_6_1[6] = {13,27,26,24,21,7};
static unsigned const b2_6_2[5] = {14,25,23,20,6};
static unsigned const b2_6_3[4] = {15,22,19,5};
static unsigned const b2_6_4[3] = {16,18,4};
static unsigned const b2_6_5[2] = {17,3};
static unsigned const b2_6_6[1] = {0};
static unsigned const* const b2_6[7] =
{b2_6_0,b2_6_1,b2_6_2,b2_6_3,b2_6_4,b2_6_5,b2_6_6};
static unsigned const b2_7_0[8] = {2,14,13,12,11,10,9,1};
static unsigned const b2_7_1[7] = {15,35,34,32,29,25,8};
static unsigned const b2_7_2[6] = {16,33,31,28,24,7};
static unsigned const b2_7_3[5] = {17,30,27,23,6};
static unsigned const b2_7_4[4] = {18,26,22,5};
static unsigned const b2_7_5[3] = {19,21,4};
static unsigned const b2_7_6[2] = {20,3};
static unsigned const b2_7_7[1] = {0};
static unsigned const* const b2_7[8] =
{b2_7_0,b2_7_1,b2_7_2,b2_7_3,b2_7_4,b2_7_5,b2_7_6,b2_7_7};
static unsigned const b2_8_0[9] = {2,16,15,14,13,12,11,10,1};
static unsigned const b2_8_1[8] = {17,44,43,41,38,34,29,9};
static unsigned const b2_8_2[7] = {18,42,40,37,33,28,8};
static unsigned const b2_8_3[6] = {19,39,36,32,27,7};
static unsigned const b2_8_4[5] = {20,35,31,26,6};
static unsigned const b2_8_5[4] = {21,30,25,5};
static unsigned const b2_8_6[3] = {22,24,4};
static unsigned const b2_8_7[2] = {23,3};
static unsigned const b2_8_8[1] = {0};
static unsigned const* const b2_8[9] =
{b2_8_0,b2_8_1,b2_8_2,b2_8_3,b2_8_4,b2_8_5,b2_8_6,b2_8_7,b2_8_8};
static unsigned const b2_9_0[10] = {2,18,17,16,15,14,13,12,11,1};
static unsigned const b2_9_1[9] = {19,54,53,51,48,44,39,33,10};
static unsigned const b2_9_2[8] = {20,52,50,47,43,38,32,9};
static unsigned const b2_9_3[7] = {21,49,46,42,37,31,8};
static unsigned const b2_9_4[6] = {22,45,41,36,30,7};
static unsigned const b2_9_5[5] = {23,40,35,29,6};
static unsigned const b2_9_6[4] = {24,34,28,5};
static unsigned const b2_9_7[3] = {25,27,4};
static unsigned const b2_9_8[2] = {26,3};
static unsigned const b2_9_9[1] = {0};
static unsigned const* const b2_9[10] =
{b2_9_0,b2_9_1,b2_9_2,b2_9_3,b2_9_4,b2_9_5,b2_9_6,b2_9_7,b2_9_8,b2_9_9};
static unsigned const b2_10_0[11] = {2,20,19,18,17,16,15,14,13,12,1};
static unsigned const b2_10_1[10] = {21,65,64,62,59,55,50,44,37,11};
static unsigned const b2_10_2[9] = {22,63,61,58,54,49,43,36,10};
static unsigned const b2_10_3[8] = {23,60,57,53,48,42,35,9};
static unsigned const b2_10_4[7] = {24,56,52,47,41,34,8};
static unsigned const b2_10_5[6] = {25,51,46,40,33,7};
static unsigned const b2_10_6[5] = {26,45,39,32,6};
static unsigned const b2_10_7[4] = {27,38,31,5};
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

int computeTriNodeIndex(int P, int i, int j)
{
  int k = P-i-j;
  if(i == P) return 0;
  if(j == P) return 1;
  if(k == P) return 2;
  if(k == 0) return 2+j; // 0-1
  if(i == 0) return 2+(P-1)+k; // 1-2
  if(j == 0) return 2+(P-1)*2+i; // 2-0
  return k*(P-1)-k*(k-1)/2+j+2*P;
}

int computeTetNodeIndex(int P, int i, int j, int k)
{
  int l = P-i-j-k;
  if(i == P) return 0;
  if(j == P) return 1;
  if(k == P) return 2;
  if(l == P) return 3;
  if(k == 0 && l == 0) return 3+j; // 0-1
  if(i == 0 && l == 0) return 3+(P-1)+k; // 1-2
  if(j == 0 && l == 0) return 3+2*(P-1)+i; // 2-0
  if(j == 0 && k == 0) return 3+3*(P-1)+l; // 0-3
  if(i == 0 && k == 0) return 3+4*(P-1)+l; // 1-3
  if(i == 0 && j == 0) return 3+5*(P-1)+l;// 2-3
  if(l == 0) return k*(P-1)-k*(k-1)/2+j+5*P-2; // 0-1-2
  if(k == 0) return l*(P-1)-l*(l-1)/2+j+5*P-2+(P-2)*(P-1)/2; // 0-1-3
  if(i == 0) return l*(P-1)-l*(l-1)/2+k+5*P-2+(P-2)*(P-1);// 1-2-3
  if(j == 0) return l*(P-1)-l*(l-1)/2+k+5*P-2+(P-2)*(P-1)*3/2; // 0-2-3
  return i-P-((i-P+1)*(i-P+2)*(i-P+3))/6+l*(P-1-i)-l*(l-1)/2+k+2*P*P+2;
}

// publically accessible access
// There is likely room for improvement, for example if a dynamic numbering
// is called on the fly, compute it, store it, and then read it off the table
int getTriNodeIndex(int P, int i, int j)
{
  // use a table if its small, otherwise dynamically generate it on the fly
  if(P <= 10)
    return crv::b2[P][i][j];
  else
    return computeTriNodeIndex(P,i,j);
}

int getTetNodeIndex(int P, int i, int j, int k)
{
  if(P <= 4)
    return crv::b3[P][i][j][k];
  else
    return computeTetNodeIndex(P,i,j,k);
}

// f is face number, see apf::tet_tri_verts or other docs
template <class T>
static void getTriFromTet(int f, int P, apf::NewArray<T>& tetNodes,
    apf::NewArray<T>& triNodes)
{
  int IJKL[4] = {0,0,0,0};
  for (int i = 0; i <= P; ++i)
    for (int j = 0; j <= P-i; ++j){
      // one of I,J,K,L = 0
      // one of I,J,K,L = P - the other two
      IJKL[apf::tet_tri_verts[f][0]] = i;
      IJKL[apf::tet_tri_verts[f][1]] = j;
      IJKL[apf::tet_tri_verts[f][2]] = P-i-j;
      triNodes[getTriNodeIndex(P,i,j)] =
        tetNodes[getTetNodeIndex(P,IJKL[0],IJKL[1],IJKL[2])];
    }
}

void getTriNodesFromTetNodes(int f, int P,
    apf::NewArray<apf::Vector3>& tetNodes,
    apf::NewArray<apf::Vector3>& triNodes)
{
  getTriFromTet(f,P,tetNodes,triNodes);
}

void getTriDetJacobianNodesFromTetDetJacobianNodes(int f, int P,
    apf::NewArray<double>& tetNodes,
    apf::NewArray<double>& triNodes)
{
  getTriFromTet(f,P,tetNodes,triNodes);
}

static unsigned const tet_tri4_f0r0[3] = {0,1,2};
static unsigned const tet_tri4_f0r1[3] = {2,0,1};
static unsigned const tet_tri4_f0r2[3] = {1,2,0};
static unsigned const tet_tri4_f1r0[3] = {2,1,0};
static unsigned const tet_tri4_f1r1[3] = {1,0,2};
static unsigned const tet_tri4_f1r2[3] = {0,2,1};

static unsigned const* const tet_tri4_f0[3] =
{tet_tri4_f0r0,tet_tri4_f0r1,tet_tri4_f0r2};
static unsigned const* const tet_tri4_f1[3] =
{tet_tri4_f1r0,tet_tri4_f1r1,tet_tri4_f1r2};
static unsigned const* const* const tet_tri4[2] = {tet_tri4_f0,tet_tri4_f1};

static unsigned const tet_tri5_f0r0[6] = {0,1,2,3,4,5};
static unsigned const tet_tri5_f0r1[6] = {5,3,0,4,1,2};
static unsigned const tet_tri5_f0r2[6] = {2,4,5,1,3,0};
static unsigned const tet_tri5_f1r0[6] = {5,4,2,3,1,0};
static unsigned const tet_tri5_f1r1[6] = {2,1,0,4,3,5};
static unsigned const tet_tri5_f1r2[6] = {0,3,5,1,4,2};

static unsigned const* const tet_tri5_f0[3] =
{tet_tri5_f0r0,tet_tri5_f0r1,tet_tri5_f0r2};
static unsigned const* const tet_tri5_f1[3] =
{tet_tri5_f1r0,tet_tri5_f1r1,tet_tri5_f1r2};
static unsigned const* const* const tet_tri5[2] = {tet_tri5_f0,tet_tri5_f1};

static unsigned const tet_tri6_f0r0[10] = {0,1,2,3,4,5,6,7,8,9};
static unsigned const tet_tri6_f0r1[10] = {9,7,4,0,8,5,1,6,2,3};
static unsigned const tet_tri6_f0r2[10] = {3,6,8,9,2,5,7,1,4,0};
static unsigned const tet_tri6_f1r0[10] = {9,8,6,3,7,5,2,4,1,0};
static unsigned const tet_tri6_f1r1[10] = {3,2,1,0,6,5,4,8,7,9};
static unsigned const tet_tri6_f1r2[10] = {0,4,7,9,1,5,8,2,6,3};

static unsigned const* const tet_tri6_f0[3] =
{tet_tri6_f0r0,tet_tri6_f0r1,tet_tri6_f0r2};
static unsigned const* const tet_tri6_f1[3] =
{tet_tri6_f1r0,tet_tri6_f1r1,tet_tri6_f1r2};
static unsigned const* const* const tet_tri6[2] = {tet_tri6_f0,tet_tri6_f1};

unsigned const* const* const* const tet_tri[7] =
{0,0,0,0,tet_tri4,tet_tri5,tet_tri6};

}
