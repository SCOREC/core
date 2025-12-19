#include "maReposition.h"
#include "maAdapt.h"
#include "maShape.h"
#include "maShapeHandler.h"
#include <apfMesh.h>
#include <apf.h>
#include <mthMatrix.h>
#include <mthAD.h>
#include <mth_def.h>
#include <apf2mth.h>

#include <sstream>
#include <cmath>
#include <functional>

namespace ma {

RepositionVertex::RepositionVertex(Adapt* a) : adapt(a), mesh(a->mesh)
{
}

void RepositionVertex::init(Entity* vertex)
{
  this->invalid.n = 0;
  this->adjEdges.n = 0;
  this->adjacentElements.setSize(0);
  this->vertex = vertex;
  this->prevPosition = getPosition(mesh, vertex);
  mesh->getAdjacent(vertex, mesh->getDimension(), adjacentElements);
  storeOldCache();
}

void RepositionVertex::storeOldCache()
{
  oldCache.clear();
  for (size_t i = 0; i < adjacentElements.getSize(); ++i) {
    if (mesh->hasTag(adjacentElements[i], adapt->qualityCache)) {
      oldCache.push_back(getCachedQuality(adapt, adjacentElements[i]));
      mesh->removeTag(adjacentElements[i], adapt->qualityCache);
    }
    else oldCache.push_back(-1);
  }
}

void RepositionVertex::findInvalid()
{
  for (size_t i = 0; i < adjacentElements.getSize(); ++i) {
    /* for now, when snapping a vertex on the boundary
    layer, ignore the quality of layer elements.
    not only do we not have metrics for this, but the
    algorithm that moves curves would need to change */
    if (getFlag(adapt, adjacentElements[i], LAYER)) continue;

    double quality;
    if (mesh->getType(adjacentElements[i]) == apf::Mesh::TET && !isTetValid(mesh, adjacentElements[i]))
      quality = -1;
    else 
      quality = adapt->shape->getQuality(adjacentElements[i]);

    if (quality < adapt->input->validQuality) invalid.e[invalid.n++] = adjacentElements[i];
  }
}

bool RepositionVertex::move(Entity* vertex, Vector target)
{
  init(vertex);
  mesh->setPoint(vertex, 0, target);
  findInvalid();
  if (invalid.n == 0) return true;
  cancel(vertex);
  return false;
}

double RepositionVertex::findWorstShape(Vector position)
{
  mesh->setPoint(vertex, 0, position);
  double worst = 1;
  for (size_t i = 0; i < adjacentElements.getSize(); ++i) {
    double quality = adapt->shape->getQuality(adjacentElements[i]);
    if (quality < worst) worst = quality;
  }
  return worst;
}

double goldenSearch(const std::function<double(Vector)> &f, Vector left, Vector right, double tol)
{
  const double M_GOLDEN_RATIO = (1.0 + std::sqrt(5.0)) / 2.0;
  double leftValue;

  do {
    Vector delta_left = left + (right - left) * M_GOLDEN_RATIO;
    Vector delta_right = right + (left - right) * M_GOLDEN_RATIO;
    leftValue = f(delta_left);

    if (leftValue < f(delta_right))
      right = delta_right;
    else
      left = delta_left;

  } while ((left-right).getLength() > tol);
  return leftValue;
}

bool RepositionVertex::moveToImproveQuality(Entity* vertex)
{
  if (mesh->getModelType(mesh->toModel(vertex)) != 3) return false; //TODO: remove limitation
  init(vertex);
  mesh->getUp(vertex, adjEdges);

  Vector center = modelCenter();
  Vector target = center + (prevPosition - center)/4;
  const auto getQuality = [this](const Vector& pos) { return this->findWorstShape(pos); };
  worstQuality = goldenSearch(getQuality, prevPosition, target, 0.000001);
  return worstQuality > adapt->input->goodQuality;
}

bool RepositionVertex::improveQuality(Entity* tet)
{
  Entity* verts[4];
  mesh->getDownward(tet, 0, verts);
  for (int i=0; i<4; i++)
    if (moveToImproveQuality(verts[i])) return true;
  return false;
}

Vector RepositionVertex::modelCenter()
{
  Model* vertModel = mesh->toModel(vertex);
  Vector avg(0,0,0);
  for (int i=0; i<adjEdges.n; i++) {
    if (vertModel != mesh->toModel(adjEdges.e[i])) continue;
    Entity* opp = getEdgeVertOppositeVert(mesh, adjEdges.e[i], vertex);
    avg += getPosition(mesh, opp);
  }
  avg = avg / adjEdges.n;
  return avg;
}

double RepositionVertex::findShortestEdge(Vector position)
{
  double quality = findWorstShape(position);
  if (quality < startingQuality && quality < adapt->input->goodQuality) return -1;

  double shortest = 1000;
  for (int i=0; i<adjEdges.n; i++) {
    double length = adapt->sizeField->measure(adjEdges.e[i]);
    if (length < shortest) shortest = length;
  }
  return shortest;
}

void RepositionVertex::moveToImproveShortEdges(Entity* vertex)
{
  int dim = mesh->getModelType(mesh->toModel(vertex));
  if (dim == 0) return;
  init(vertex);
  
  mesh->getUp(vertex, adjEdges);
  this->startingQuality = findWorstShape(prevPosition);
  Vector center = modelCenter();
  const auto getLength = [this](const Vector& pos) { return this->findShortestEdge(pos); };
  goldenSearch(getLength, prevPosition, center, 0.000001);
}

void RepositionVertex::cancel(Entity* vertex)
{
  PCU_ALWAYS_ASSERT(this->vertex == vertex);
  mesh->setPoint(vertex, 0, prevPosition);
}

apf::Up& RepositionVertex::getInvalid()
{
  return invalid;
}

typedef mth::AD<double, 3> AD;
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
