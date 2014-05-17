/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#include <cfloat>
#include "maMesh.h"
#include "maSize.h"
#include "maAdapt.h"
#include "maShapeHandler.h"

namespace ma {

double measureTriQuality(Mesh* m, SizeField* f, Entity* tri)
{
  Entity* e[3];
  m->getDownward(tri,1,e);
  double l[3];
  for (int i=0; i < 3; ++i)
    l[i] = f->measure(e[i]);
  double A = f->measure(tri);
  double s = 0;
  for (int i=0; i < 3; ++i)
    s += l[i]*l[i];
  return 48*(A*A)/(s*s);
}

/* applies the mean ratio cubed formula from Li's thesis */
double measureTetQuality(Mesh* m, SizeField* f, Entity* tet)
{
  Entity* e[6];
  m->getDownward(tet,1,e);
  double l[6];
  for (int i=0; i < 6; ++i)
    l[i] = f->measure(e[i]);
  double V = f->measure(tet);
  double s=0;
  for (int i=0; i < 6; ++i)
    s += l[i]*l[i];
  if (V < 0)
    return -15552*(V*V)/(s*s*s);
  return 15552*(V*V)/(s*s*s);
}

double measureElementQuality(Mesh* m, SizeField* f, Entity* e)
{
  typedef double (*MeasureQualityFunction)(Mesh*,SizeField*,Entity*);
  static MeasureQualityFunction table[TYPES] =
  {0
  ,0
  ,measureTriQuality
  ,0
  ,measureTetQuality
  ,0
  ,0
  ,0};
  return table[m->getType(e)](m,f,e);
}

double getWorstQuality(Adapt* a, Entity** e, size_t n)
{
  Mesh* m = a->mesh;
  ShapeHandler* sh = a->shape;
  double worst = 1;
  for (size_t i=0; i < n; ++i)
  {
    double quality = sh->getQuality(e[i]);
    if (quality < worst) worst = quality;
  }
  return worst;
}

double getWorstQuality(Adapt* a, EntityArray& e)
{
  return getWorstQuality(a,&(e[0]),e.getSize());
}

/* applies the same measure as measureTetQuality
   but works directly off the points. */
double measureLinearTetQuality(Vector xyz[4])
{
  Matrix J;
  J[0] = xyz[1] - xyz[0];
  J[1] = xyz[2] - xyz[0];
  J[2] = xyz[3] - xyz[0];
  double j = apf::det(J);
  double V = j/6;
  double s = 0;
  for (int i = 0; i < 6; ++i) {
    int const* ev = apf::tet_edge_verts[i];
    double l = (xyz[ev[1]] - xyz[ev[0]]).getLength();
    s += l * l;
  }
  if (V < 0)
    return -15552*(V*V)/(s*s*s);
  return 15552*(V*V)/(s*s*s);
}

/* helper for measureBezierTetQuality only.
   hardcoded the only inputs for speed.*/
static int factorial(int num)
{
  static int const table[4] = {1,1,2,6};
  return table[num];
}

/* code provided by Quikai Lu for computing
   the shape quality of a 2nd-order Bezier tetrahedron (10-point).
   The output value is in (-inf,0] for invalid elements
   and in (0,1] for valid elements, with values closer to 1
   being better quality.
   Note that this is a curving quality measure only, it outputs 1 for
   all straight-sided elements. */
double measureBezierTetQuality(Vector xyz[10])
{

  bool isValid = true;
  /* if the ST part is valid, further compute the curved (CR) part */

  Vector tetCtrlPts[3][3][3];

  Vector du[2][2][2];
  Vector dv[2][2][2];
  Vector dw[2][2][2];

  double jac[4][4][4];

  int i, j, k, l, m, n, r, s, t, d1, d2, d3, d4, d5, d6;
  static int const degree = 2;
  double minDetJ = DBL_MAX;
  double maxDetJ = -DBL_MAX;

  // vertices
  tetCtrlPts[0][0][0] = xyz[0];
  tetCtrlPts[2][0][0] = xyz[1];
  tetCtrlPts[0][2][0] = xyz[2];
  tetCtrlPts[0][0][2] = xyz[3];

  // edges
  tetCtrlPts[1][0][0] = xyz[4];
  tetCtrlPts[1][1][0] = xyz[5];
  tetCtrlPts[0][1][0] = xyz[6];
  tetCtrlPts[0][0][1] = xyz[7];
  tetCtrlPts[1][0][1] = xyz[8];
  tetCtrlPts[0][1][1] = xyz[9];

  //calculate the derivatives
  for(i = 0; i < degree; i++) {
    for(j = 0; j < (degree - i); j++) {
      for(k = 0; k < (degree - i - j); k++) {
        du[i][j][k] = tetCtrlPts[i+1][j][k] - tetCtrlPts[i][j][k];
        dv[i][j][k] = tetCtrlPts[i][j+1][k] - tetCtrlPts[i][j][k];
        dw[i][j][k] = tetCtrlPts[i][j][k+1] - tetCtrlPts[i][j][k];
      }
    }
  }

  //clear the jacobian 
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4-i; j++) {
      for (k = 0; k < 4-i-j; k++) { 
        jac[i][j][k] = 0.0;
      }
    }
  }

  //calculate the box products
  static double const degree2 = degree * degree * degree;
  double detj, coeff;

  /* we await a hero to slay this nested monster */
  for (i = 0; i < degree; i++) {
    d1 = degree - i;
    for (j = 0; j < d1; j++) {
      d2 = d1 - j;
      for(k = 0; k < d2; k++){
        for (l = 0; l < degree; l++) {
          d3 = degree - l;
          for (m = 0; m < d3; m++) {
            d4 = d3 - m;
            for(n = 0; n < d4; n++){
              for(r = 0; r < degree; r++){
                d5 = degree - r;
                for(s = 0; s < d5; s++){
                  d6 = d5 - s;
                  for(t = 0; t< d6; t++){
                    detj = apf::cross(du[i][j][k], dv[l][m][n]) * dw[r][s][t];
                    jac[i+l+r][j+m+s][k+n+t] += detj * degree2;
                  }
                }
              }
            }
          }
        }
      }
    }
  } 

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4-i; j++) {
      for (k = 0; k < 4-i-j; k++) { 
        coeff =  (double)factorial(3*(degree-1))
                /(double)(factorial(i)*
                          factorial(j)*
                          factorial(k)*
                          factorial(3*(degree-1)-i-j-k));
        jac[i][j][k] /= coeff;
        if(jac[i][j][k] <= 0.0)
          isValid = false;
        minDetJ = std::min(minDetJ, jac[i][j][k]);
        maxDetJ = std::max(maxDetJ, jac[i][j][k]);
      }
    }
  }

  if(isValid)
    return minDetJ / maxDetJ;
  return minDetJ;
}

static Vector getBezierEdgePoint(
    Vector const& l1,
    Vector const& l2,
    Vector const& l3)
{
  return ((l2 * 4.0) - (l1 + l3)) / 2.0;
}

static void convertQuadraticTetToBezier(Vector xyz[10])
{
  Vector* ep = xyz + 4;
  for (int i = 0; i < 6; ++i) {
    int const* ev = apf::tet_edge_verts[i];
    ep[i] = getBezierEdgePoint(xyz[ev[0]],ep[i],xyz[ev[1]]);
  }
}

double measureQuadraticTetQuality(Vector xyz[10])
{
  double lq = measureLinearTetQuality(xyz);
  if (lq <= 0)
    return lq;
  convertQuadraticTetToBezier(xyz);
  double qq = measureBezierTetQuality(xyz);
  if (qq <= 0)
    return qq;
  /* since we don't really have correction procedures
     for the curved edges, lets just use the linear
     quality measure for valid tets, and employ
     the bezier code to ensure that the curved version
     is valid at all points */
  return lq;
}

/* serious theory warning: this function does NOT
   measure quality in metric space, it can only
   do it in real space.
   This means it is only good for isotropic size
   fields right now, where the size field does
   not affect quality */
double measureQuadraticTetQuality(Mesh* m, Entity* tet)
{
  Vector xyz[10];
  Vector* ep = xyz + 4;
  Entity* v[4];
  m->getDownward(tet,0,v);
  for (int i = 0; i < 4; ++i)
    m->getPoint(v[i],0,xyz[i]);
  Entity* e[6];
  m->getDownward(tet,1,e);
  for (int i = 0; i < 6; ++i)
    m->getPoint(e[i],0,ep[i]);
  return measureQuadraticTetQuality(xyz);
}

}

