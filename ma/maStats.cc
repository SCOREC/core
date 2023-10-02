/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "maStats.h"
#include "maAdapt.h"

namespace ma {

void getLinearQualitiesInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &linearQualities)
{
  ma::Entity* e;
  ma::Iterator* it;
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    if (! m->isOwned(e))
      continue;
    if (! apf::isSimplex(m->getType(e))) // ignore non-simplex elements
      continue;
    double lq = ma::measureElementQuality(m, sf, e);
    if (m->getDimension() == 2)
      lq = (lq > 0) ? std::sqrt(lq) : -std::sqrt(-lq);
    else
      lq = cbrt(lq);
    linearQualities.push_back(lq);
  }
  m->end(it);
}

void getEdgeLengthsInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths)
{
  ma::Entity* e;
  ma::Iterator* it;
  it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    if (! m->isOwned(e))
      continue;
    edgeLengths.push_back(sf->measure(e));
  }
  m->end(it);
}

void getLinearQualitiesInPhysicalSpace(ma::Mesh* m,
    std::vector<double> &linearQualities)
{
  ma::Entity* e;
  ma::Iterator* it;
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    double lq;
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
      if (A < 0)
	lq = -48. * (A*A) / (s*s);
      else
	lq =  48. * (A*A) / (s*s);
    }
    else if (m->getType(e) == apf::Mesh::TET){
      ma::Vector p[4];
      ma::getVertPoints(m, e, p);
      lq = ma::measureLinearTetQuality(p);
    }
    if (m->getDimension() == 2)
      lq = (lq > 0) ? std::sqrt(lq) : -std::sqrt(-lq);
    else
      lq = cbrt(lq);
    linearQualities.push_back(lq);
  }
  m->end(it);
}

void getEdgeLengthsInPhysicalSpace(ma::Mesh* m,
    std::vector<double> &edgeLengths)
{
  ma::Entity* e;
  ma::Iterator* it;
  it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    SizeField* sf = new IdentitySizeField(m);
    edgeLengths.push_back(sf->measure(e));
  }
  m->end(it);
}

void getStatsInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities)
{
  // linear qualities
  getLinearQualitiesInMetricSpace(m, sf, linearQualities);
  // edge lengths
  getEdgeLengthsInMetricSpace(m, sf, edgeLengths);
}

void getStatsInPhysicalSpace(ma::Mesh* m,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities)
{
  // linear qualities
  getLinearQualitiesInPhysicalSpace(m, linearQualities);
  // edge lengths
  getEdgeLengthsInPhysicalSpace(m, edgeLengths);
}

/** \brief Measures mesh statistics given a sizeField
  \details  quantities measured are:
  1) normalized edge length
  2) linear quality
  The values can be computed in both metric (if inMetric = true)
  and physical (if inMetric = false) spaces.*/
void stats(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    bool inMetric)
{
  edgeLengths.clear();
  linearQualities.clear();

  if (inMetric)
    getStatsInMetricSpace(m, sf, edgeLengths, linearQualities);
  else
    getStatsInPhysicalSpace(m, edgeLengths, linearQualities);
}

}
