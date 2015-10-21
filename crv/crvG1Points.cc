/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvBezier.h"
#include "crvTables.h"
#include <cassert>

namespace crv {

static void getGregoryTriangleTransform(int P, apf::NewArray<double> & c)
{
  apf::NewArray<double> d;
  getTransformationCoefficients(P,apf::Mesh::TRIANGLE,d);

  int nbBezier = (P-1)*(P-2)/2;
  int niBezier = (P+1)*(P+2)/2;

  int nb = 6;
  int ni = 6+3*P;
  c.allocate(ni*nb);

  int map[3] = {1,2,0};
  // copy the bezier point locations
  for(int i = 0; i < nbBezier; ++i){
    for(int j = 0; j < niBezier; ++j)
      c[i*ni+j] = d[i*niBezier+j];
    for(int j = niBezier; j < ni; ++j)
      c[i*ni+j] = 0.;
  }
  if(P == 3){
    for(int i = nbBezier; i < nb; ++i){
      for(int j = 0; j < niBezier; ++j)
        c[i*ni+j] = d[j];
      for(int j = niBezier; j < ni; ++j)
        c[i*ni+j] = 0.;
    }
  }
  if(P == 4){
    for(int i = nbBezier; i < nb; ++i){
      for(int j = 0; j < niBezier; ++j)
        c[i*ni+j] = d[map[i-nbBezier]*niBezier+j];
      for(int j = niBezier; j < ni; ++j)
        c[i*ni+j] = 0.;
    }
  }
}

static void getGregoryTetTransform(int P, apf::NewArray<double> & c)
{
  assert(P == 4 && getBlendingOrder(apf::Mesh::TET) == 0);
  double t4[47] = {
      -0.665492638178598,-0.665492638178598,-0.665492638178598,-0.665492638178598,
      0.697909481209196,0.496340368840329,0.697909481209197,0.697909481209196,
      0.496340368840329,0.697909481209196,0.697909481209196,0.49634036884033,
      0.697909481209196,0.697909481209196,0.496340368840329,0.697909481209196,
      0.697909481209196,0.496340368840329,0.697909481209196,0.697909481209196,
      0.496340368840329,0.697909481209196,
      -0.764902170896025,-0.764902170896025,-0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,-0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,-0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,-0.764902170896025,-0.764902170896025,
      -0.764902170896025,-0.764902170896025,
      10.6666666666667,
  };

  c.allocate(47);
  for (int j = 0; j < 47; ++j)
    c[j] = t4[j];
}

static void getBlendedGregoryTriangleTransform(int P, int blend,
    apf::NewArray<double> & c)
{
  apf::NewArray<double> d;
  getBlendedTransformationCoefficients(P,blend,apf::Mesh::TRIANGLE,d);

  int nbBezier = (P-1)*(P-2)/2;
  int niBezier = (P+1)*(P+2)/2-nbBezier;

  int nb = 6;
  int ni = 3*P;
  c.allocate(ni*nb);

  int map[3] = {1,2,0};
  // copy the bezier point locations
  for(int i = 0; i < nbBezier; ++i){
    for(int j = 0; j < niBezier; ++j)
      c[i*ni+j] = d[i*niBezier+j];
    for(int j = niBezier; j < ni; ++j)
      c[i*ni+j] = 0.;
  }
  if(P == 3){
    for(int i = nbBezier; i < nb; ++i){
      for(int j = 0; j < niBezier; ++j)
        c[i*ni+j] = d[j];
      for(int j = niBezier; j < ni; ++j)
        c[i*ni+j] = 0.;
    }
  }
  if(P == 4){
    for(int i = nbBezier; i < nb; ++i){
      for(int j = 0; j < niBezier; ++j)
        c[i*ni+j] = d[map[i-nbBezier]*niBezier+j];
      for(int j = niBezier; j < ni; ++j)
        c[i*ni+j] = 0.;
    }
  }
}

static void getBlendedGregoryTetTransform(int P, int blend,
    apf::NewArray<double> & c)
{
  assert(P == 4 && getBlendingOrder(apf::Mesh::TET) == 0);
  double t4_1[46] = {
      1.921296296296296,1.921296296296296,1.921296296296296,1.921296296296296,
      -0.7098765432098767,-1.064814814814814,-0.7098765432098767,-0.7098765432098767,
      -1.064814814814814,-0.7098765432098767,-0.7098765432098767,-1.064814814814814,
      -0.7098765432098767,-0.7098765432098767,-1.064814814814814,-0.7098765432098767,
      -0.7098765432098767,-1.064814814814814,-0.7098765432098767,-0.7098765432098767,
      -1.064814814814814,-0.7098765432098767,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593,
      0.342592592592593,0.342592592592593,0.342592592592593};
  double t4_2[46] = {
      0.3472222222222222,0.3472222222222222,0.3472222222222222,0.3472222222222222,
      -0.2407407407407403,-0.3611111111111111,-0.2407407407407404,-0.2407407407407404,
      -0.3611111111111111,-0.2407407407407404,-0.2407407407407404,-0.3611111111111111,
      -0.2407407407407404,-0.2407407407407404,-0.3611111111111111,-0.2407407407407404,
      -0.2407407407407404,-0.3611111111111111,-0.2407407407407404,-0.2407407407407404,
      -0.3611111111111111,-0.2407407407407404,
      0.388888888888888,0.3888888888888888,0.3888888888888888,0.3888888888888888,
      0.3888888888888888,0.3888888888888888,
      0.3888888888888888,0.3888888888888888,0.3888888888888888,0.3888888888888888,
      0.3888888888888888,0.3888888888888888,
      0.388888888888888,0.3888888888888888,0.3888888888888888,0.3888888888888888,
      0.3888888888888888,0.3888888888888888,
      0.3888888888888888,0.3888888888888888,0.3888888888888888,0.3888888888888888,
      0.3888888888888888,0.3888888888888888};
  double* t4[2] = {t4_1,t4_2};
  c.allocate(46);
  for (int j = 0; j < 46; ++j)
    c[j] = t4[blend-1][j];
}

void getGregoryTransformationCoefficients(int P, int type,
    apf::NewArray<double>& c){
  assert(P == 3 || P == 4);
  if(type == apf::Mesh::EDGE)
    getTransformationCoefficients(P,apf::Mesh::EDGE,c);
  else if(type == apf::Mesh::TRIANGLE)
    getGregoryTriangleTransform(P,c);
  else if(type == apf::Mesh::TET)
    getGregoryTetTransform(P,c);
}

void getGregoryBlendedTransformationCoefficients(int P, int blend, int type,
    apf::NewArray<double>& c){

  if(type == apf::Mesh::TRIANGLE)
    getBlendedGregoryTriangleTransform(P,blend,c);
  else if(type == apf::Mesh::TET)
    getBlendedGregoryTetTransform(P,blend,c);
}

}
