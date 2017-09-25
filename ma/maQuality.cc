/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#include <cfloat>
#include <pcu_util.h>
#include <cstdlib>
#include "maMesh.h"
#include "maSize.h"
#include "maAdapt.h"
#include "maShapeHandler.h"
#include "maShape.h"
#include <apfGeometry.h>

namespace ma {

bool areTetsValid(Mesh* m, EntityArray& tets)
{
  Vector v[4];
  size_t n = tets.getSize();
  for (size_t i = 0; i < n; ++i){
    ma::getVertPoints(m,tets[i],v);
    if ((cross((v[1] - v[0]), (v[2] - v[0])) * (v[3] - v[0])) < 0)
      return false;
  }
  return true;
}

class FixedMetricIntegrator : public apf::Integrator
{
  public:
    FixedMetricIntegrator(Mesh* inMesh, const Matrix& inQ):
      Integrator(1),
      measurement(0),
      mesh(inMesh),
      Q(inQ)
    {
      dimension = 0;
      elem = 0;
    }
    virtual void inElement(apf::MeshElement* me)
    {
      dimension = apf::getDimension(me);
      elem = apf::createElement(mesh->getCoordinateField(), me);
    }
    virtual void outElement()
    {
      apf::destroyElement(elem);
    }
    virtual void atPoint(Vector const& p, double w, double)
    {
      Matrix J;
      apf::getJacobian(apf::getMeshElement(elem),p,J);
/* transforms the rows of J, the differential tangent vectors,
   into the metric space, then uses the generalized determinant */
      double dV2 = apf::getJacobianDeterminant(J*Q,dimension);
      measurement += w*dV2;
    }
    double measurement;
  private:
    Mesh* mesh;
    Matrix Q;
    int dimension;
    apf::Element* elem;
};


double qMeasure(Mesh* mesh, Entity* e, const Matrix& Q)
{
  FixedMetricIntegrator integrator(mesh, Q);
  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  integrator.process(me);
  apf::destroyMeshElement(me);
  return integrator.measurement;
}

static Matrix getMetricWithMaxJacobean(Mesh* m, SizeField* sf,
    Entity* e)
{
  int dim = m->getDimension();
  int type = m->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TRIANGLE ||
		    type == apf::Mesh::TET);
  Downward dv;
  int nd = m->getDownward(e, 0, dv);

  Matrix Q;
  double maxJ = -1.0;

  for (int i = 0; i < nd; i++) {
    apf::MeshElement* me = createMeshElement(m, dv[i]);
    Matrix currentQ;
    sf->getTransform(me, Vector(0.0, 0.0, 0.0), currentQ);
    double currentJ = apf::getJacobianDeterminant(currentQ, dim);
    if (currentJ > maxJ) {
      maxJ = currentJ;
      Q = currentQ;
    }
    apf::destroyMeshElement(me);
  }
  return Q;
}

double measureTriQuality(Mesh* m, SizeField* f, Entity* tri, bool useMax)
{
  /* By default, we are using Q at the center of the tri.
   * If useMax is true metric at a (downward) vertex with the
   * largest determinant is used.
   * Note: In the future we may want to used average of Q over the tri */
  Matrix Q;
  if (useMax)
    Q = getMetricWithMaxJacobean(m, f, tri);
  else {
    apf::MeshElement* me = createMeshElement(m, tri);
    Vector xi(1./3., 1./3., 1./3.);
    f->getTransform(me, xi, Q);
    apf::destroyMeshElement(me);
  }

  Entity* e[3];
  m->getDownward(tri,1,e);
  double l[3];
  for (int i=0; i < 3; ++i)
    l[i] = qMeasure(m, e[i], Q);
  double A = qMeasure(m, tri, Q);
  double s = 0;
  for (int i=0; i < 3; ++i)
    s += l[i]*l[i];
  return 48*(A*A)/(s*s);
}

/* applies the mean ratio cubed formula from Li's thesis */
double measureTetQuality(Mesh* m, SizeField* f, Entity* tet, bool useMax)
{
  /* By default, we are using Q at the center of the tet.
   * If useMax is true metric at a (downward) vertex with the
   * largest determinant is used.
   * Note: In the future we may want to used average of Q over the tet */
  Matrix Q;
  if (useMax)
    Q = getMetricWithMaxJacobean(m, f, tet);
  else {
    apf::MeshElement* me = createMeshElement(m, tet);
    Vector xi(0.25, 0.25, 0.25);
    f->getTransform(me, xi, Q);
    apf::destroyMeshElement(me);
  }

  Entity* e[6];
  m->getDownward(tet,1,e);
  double l[6];
  for (int i=0; i < 6; ++i)
    l[i] = qMeasure(m, e[i], Q);
  double V = qMeasure(m, tet, Q);
  double s=0;
  for (int i=0; i < 6; ++i)
    s += l[i]*l[i];
  if (V < 0)
    return -15552*(V*V)/(s*s*s);
  return 15552*(V*V)/(s*s*s);
}

double measureElementQuality(Mesh* m, SizeField* f, Entity* e, bool useMax)
{
  typedef double (*MeasureQualityFunction)(Mesh*,SizeField*,Entity*,bool);
  static MeasureQualityFunction table[apf::Mesh::TYPES] =
  {0
  ,0
  ,measureTriQuality
  ,0
  ,measureTetQuality
  ,0
  ,0
  ,0};
  return table[m->getType(e)](m,f,e,useMax);
}

double getWorstQuality(Adapt* a, Entity** e, size_t n)
{
  PCU_ALWAYS_ASSERT(n);
  ShapeHandler* sh = a->shape;
  double worst = sh->getQuality(e[0]);
  for (size_t i = 1; i < n; ++i) {
    double quality = sh->getQuality(e[i]);
    if (quality < worst)
      worst = quality;
  }
  return worst;
}

double getWorstQuality(Adapt* a, EntityArray& e)
{
  return getWorstQuality(a, &(e[0]), e.getSize());
}

bool hasWorseQuality(Adapt* a, EntityArray& e, double qualityToBeat)
{
  size_t n = e.getSize();
  ShapeHandler* sh = a->shape;
  for (size_t i = 0; i < n; ++i) {
    double quality = sh->getQuality(e[i]);
    if (quality < qualityToBeat)
      return true;
  }
  return false;
}

/* applies the same measure as measureTetQuality
   but works directly off the points. */
double measureLinearTetQuality(Vector xyz[4])
{
  Matrix J;
  J[0] = xyz[1] - xyz[0];
  J[1] = xyz[2] - xyz[0];
  J[2] = xyz[3] - xyz[0];
  double j = apf::getDeterminant(J);
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

static int unrotate_prism_diagonal_code(int code, int rot)
{
  static int const shift_table[6] = {0,1,2,2,0,1};
  int out = 0;
  for (int i = 0; i < 3; ++i)
    if (code & (1 << i))
      out |= (1 << ((i + shift_table[rot]) % 3));
  return out;
}

bool isPrismOk(apf::Mesh* m, Entity* e,
    int* good_diagonal_codes)
{
  Entity* v[6];
  m->getDownward(e, 0, v);
  Vector p[6];
  for (int i = 0; i < 6; ++i)
    m->getPoint(v[i], 0, p[i]);
  bool all_good = true;
  if (good_diagonal_codes)
    *good_diagonal_codes = 0xFF;
  for (int i = 0; i < 6; ++i) {
    int const* new_to_old = prism_rotation[i];
    apf::Plane pl = apf::Plane::fromPoints(
        p[new_to_old[0]],
        p[new_to_old[1]],
        p[new_to_old[5]]);
    if (pl.distance(p[new_to_old[3]]) <= 0) {
      all_good = false;
      if (good_diagonal_codes) {
        int bad_code = unrotate_prism_diagonal_code(5, i);
        *good_diagonal_codes &= ~(1 << bad_code);
      }
    }
    if (pl.distance(p[new_to_old[4]]) <= 0) {
      all_good = false;
      if (good_diagonal_codes) {
        int bad_code = unrotate_prism_diagonal_code(4, i);
        *good_diagonal_codes &= ~(1 << bad_code);
      }
    }
    if (pl.distance(p[new_to_old[2]]) >= 0) {
      all_good = false;
      if (good_diagonal_codes) {
        int bad_code = unrotate_prism_diagonal_code(5, i);
        *good_diagonal_codes &= ~(1 << bad_code);
        bad_code = unrotate_prism_diagonal_code(4, i);
        *good_diagonal_codes &= ~(1 << bad_code);
      }
    }
  }
  return all_good;
}

bool isPyramidOk(apf::Mesh* m, Entity* e,
    int* good_rotation)
{
  Entity* v[5];
  m->getDownward(e, 0, v);
  Vector p[5];
  for (int i = 0; i < 5; ++i)
    m->getPoint(v[i], 0, p[i]);
  if (good_rotation)
    *good_rotation = -1;
  bool all_good = true;
  for (int i = 0; i < 2; ++i) {
    int const* new_to_old = pyramid_rotation[i];
    apf::Plane pl = apf::Plane::fromPoints(
        p[new_to_old[0]],
        p[new_to_old[2]],
        p[new_to_old[4]]);
    if (pl.distance(p[new_to_old[1]]) <= 0) {
      all_good = false;
      continue;
    }
    if (pl.distance(p[new_to_old[3]]) >= 0) {
      all_good = false;
      continue;
    }
    if (good_rotation)
      *good_rotation = i;
  }
  return all_good;
}

bool isLayerElementOk(Mesh* m, Entity* e)
{
  int type = m->getType(e);
  if (type == apf::Mesh::PYRAMID)
    return isPyramidOk(m, e);
  if (type == apf::Mesh::PRISM)
    return isPrismOk(m, e);
  abort();
  return false;
}

double getInsphere(Mesh* m, Entity* e)
{
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);

  // Insphere r of a tet computed by the forumla at
  // http://maths.ac-noumea.nc/polyhedr/stuff/tetra_sf_.htm
  // a, b, c, d are the four points of the tet
  // N_abc = (b-a) x (c-a) (x is the cross product)
  //              a_x a_y a_z 1
  // alpha = det( b_x b_y b_z 1 )
  //              a_x a_y a_z 1
  //              a_x a_y a_z 1
  // r = |alpha| / (||N_abc||+||N_abd||+||N_acd||+||N_bcd||)
  Entity* v[4];
  m->getDownward(e, 0, v);

  apf::Matrix<4,3> x;
  apf::Matrix<4,4> matrix;

  for (int i=0; i < 4; ++i) {
    x[i] = getPosition(m, v[i]);
    for (int j = 0; j < 3; ++j)
      matrix[i][j] = x[i][j];
    matrix[i][3] = 1;
  }

  double l = apf::cross(x[1]-x[0], x[2]-x[0]).getLength()
           + apf::cross(x[1]-x[0], x[3]-x[0]).getLength()
           + apf::cross(x[2]-x[0], x[3]-x[0]).getLength()
           + apf::cross(x[2]-x[1], x[3]-x[1]).getLength();

  return std::abs(apf::getDeterminant(matrix)) / l;
}

}

