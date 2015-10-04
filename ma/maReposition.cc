#include "maReposition.h"
#include <apfMesh.h>
#include <apf.h>
#include <mthMatrix.h>
#include <mthAD.h>
#include <mth_def.h>
#include <apf2mth.h>

#include <sstream>
#include <cmath>

namespace ma {

typedef mth::AD<3> AD;
typedef mth::Vector<AD, 3> ADVec;

static AD tet_volume(ADVec v[4])
{
  return (mth::cross((v[1] - v[0]), (v[2] - v[0])) * (v[3] - v[0])) / AD(6);
}

static AD triangle_area(ADVec v[3])
{
  return mth::norm(mth::cross((v[1] - v[0]), (v[2] - v[0]))) / AD(2.);
}

#define PERFECT_TRIANGLE_QUALITY (std::sqrt(3.0) / 4.0)
#define CUBE(x) ((x)*(x)*(x))
#define PERFECT_TET_QUALITY \
  ((std::sqrt(2.0) / 12.0) / CUBE(std::sqrt(PERFECT_TRIANGLE_QUALITY)))

static AD tet_quality(ADVec v[4])
{
  AD sum_asq(0);
  for (int i = 0; i < 4; ++i) {
    ADVec tri_v[3];
    for (int j = 0; j < 3; ++j)
      tri_v[j] = v[apf::tet_tri_verts[i][j]];
    AD a = triangle_area(tri_v);
    sum_asq = sum_asq + (a * a);
  }
  AD arms = mth::sqrt(sum_asq / 4.);
  AD vol = tet_volume(v);
  AD root_arms = mth::sqrt(arms);
  AD quality = vol / CUBE(root_arms);
  return quality / PERFECT_TET_QUALITY;
}

static AD tet_entity_quality(Mesh* m, Entity* tet, Entity* v)
{
  Entity* tv[4];
  m->getDownward(tet, 0, tv);
  ADVec tx[4];
  for (int i = 0; i < 4; ++i) {
    ma::Vector vx;
    m->getPoint(tv[i], 0, vx);
    for (unsigned j = 0; j < 3; ++j)
      tx[i][j] = AD(vx[j]);
    if (tv[i] == v)
      for (unsigned j = 0; j < 3; ++j)
        tx[i][j].diff(j);
  }
  return tet_quality(tx);
}

static AD min_cavity_quality(Mesh* m, apf::Adjacent& tets, Entity* v)
{
  AD min_qual(1.);
  for (size_t i = 0; i < tets.getSize(); ++i) {
    AD tet_qual = tet_entity_quality(m, tets[i], v);
    if (tet_qual.val() < min_qual.val())
      min_qual = tet_qual;
  }
  return min_qual;
}

bool repositionVertex(Mesh* m, Entity* v,
    int max_iters, double initial_speed)
{
  apf::Adjacent tets;
  m->getAdjacent(v, 3, tets);
  double speed = initial_speed;
  AD min_qual = min_cavity_quality(m, tets, v);
  double prev_qual = min_qual.val();
  for (int iter = 0; iter < max_iters; ++iter) {
    prev_qual = min_qual.val();
    double step = (1.0 - min_qual.val()) * speed;
    Vector grad(
        min_qual.dx(0),
        min_qual.dx(1),
        min_qual.dx(2));
    double grad_mag = grad * grad;
    if (!grad_mag)
      break;
    Vector motion = grad * (step / grad_mag);
    Vector vx;
    m->getPoint(v, 0, vx);
    vx += motion;
    m->setPoint(v, 0, vx);
    min_qual = min_cavity_quality(m, tets, v);
    /* diverging, reduce speed */
    if (min_qual.val() < prev_qual)
      speed /= 2;
  }
  return min_qual.val() > 0;
}

}
