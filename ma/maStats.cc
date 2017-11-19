/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "maStats.h"

namespace ma {

static bool isSimplexMesh(ma::Mesh* m)
{
  int count = 0;
  if (m->getDimension() == 2)
    count = apf::countEntitiesOfType(m, apf::Mesh2::QUAD);
  else
    count = apf::countEntitiesOfType(m, apf::Mesh2::HEX) +
            apf::countEntitiesOfType(m, apf::Mesh2::PRISM) +
            apf::countEntitiesOfType(m, apf::Mesh2::PYRAMID);
  return (count == 0);
}

void getStatsInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities)
{
  PCU_ALWAYS_ASSERT_VERBOSE(isSimplexMesh(m),
      "expecting an all simplex mesh!");

  ma::Entity* e;
  ma::Iterator* it;

  // linear qualities
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    if (! m->isOwned(e))
      continue;
    linearQualities.push_back(ma::measureElementQuality(m, sf, e));
  }
  m->end(it);

  // edge lengths
  it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    if (! m->isOwned(e))
      continue;
    apf::MeshElement* me = createMeshElement(m, e);
    ma::Matrix Q;
    sf->getTransform(me, ma::Vector(0., 0., 0.), Q);
    apf::destroyMeshElement(me);
    edgeLengths.push_back(ma::qMeasure(m, e, Q));
  }
  m->end(it);
}

void getStatsInPhysicalSpace(ma::Mesh* m,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities)
{
  PCU_ALWAYS_ASSERT_VERBOSE(isSimplexMesh(m),
      "expecting an all simplex mesh!");
  ma::Entity* e;
  ma::Iterator* it;

  // linear qualities
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    if (m->getType(e) == apf::Mesh::TRIANGLE) {
      ma::Vector p[3];
      ma::getVertPoints(m, e, p);
      double l[3];
      for (int i = 0; i < 3; i++)
	l[i] = (p[(i+1)%3] - p[i]).getLength();
      double A = 0.5 * apf::cross(p[1] - p[0], p[2] - p[0])[2];
      double s = 0;
      for (int i = 0; i < 3; i++)
	s += l[i] * l[i];
      double lq;
      if (A < 0)
	lq = -48. * (A*A) / (s*s);
      else
	lq =  48. * (A*A) / (s*s);
      linearQualities.push_back(lq);
    }
    else if (m->getType(e) == apf::Mesh::TET){
      ma::Vector p[4];
      ma::getVertPoints(m, e, p);
      linearQualities.push_back(ma::measureLinearTetQuality(p));
    }
  }
  m->end(it);

  // edge lengths
  it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    apf::MeshElement* me = createMeshElement(m, e);
    ma::Matrix Q = ma::Matrix(1.0, 0.0, 0.0,
			      0.0, 1.0, 0.0,
			      0.0, 0.0, 1.0);;
    apf::destroyMeshElement(me);
    edgeLengths.push_back(ma::qMeasure(m, e, Q));
  }
  m->end(it);
}

/** \brief Measures mesh statistics using and adapt input object
  \details  quantities measured are:
  1) normalized edge length
  2) linear quality
  The values can be computed in both metric (if inMetric = true)
  and physical (if inMetric = false) spaces.*/
void stats(ma::Input* in,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    bool inMetric)
{
  PCU_ALWAYS_ASSERT_VERBOSE(isSimplexMesh(in->mesh),
      "ma::stats expects and all simplex mesh!");
  edgeLengths.clear();
  linearQualities.clear();

  ma::validateInput(in);
  Adapt* a = new Adapt(in);

  ma::Mesh* m = a->mesh;
  ma::SizeField* sf = a->sizeField;

  if (inMetric)
    getStatsInMetricSpace(m, sf, edgeLengths, linearQualities);
  else
    getStatsInPhysicalSpace(m, edgeLengths, linearQualities);

  delete a;
  delete in;
}

}
