/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "maStats.h"
#include "maAdapt.h"
#include <numeric>
#include <stdexcept>
#include <algorithm>

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
    double lq = 0;
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

std::vector<int> printHistogramData(std::string name, std::vector<double> input, double min, double max, Mesh* m)
{
  const int nbins = 10;
  std::vector<int> count(nbins, 0);
  const double bin_size = (max-min)/(nbins*1.0);
  double inputMax = 0;
  double inputMin = max;

  for (size_t i = 0; i < input.size(); ++i) {
    if (std::isnan(input[i])) continue;
    if (input[i] > inputMax)
      inputMax = input[i];
    if (input[i] < inputMin)
      inputMin = input[i];
    int bin = (int)std::floor((input[i] - min)/bin_size);
    if (bin >= nbins) bin = nbins - 1;
    if (bin < 0) bin = 0;
    count[bin] += 1;
  }

  inputMin = m->getPCU()->Min<double>(inputMin);
  inputMax = m->getPCU()->Max<double>(inputMax);
  for (int i = 0; i < nbins; ++i) count[i] = m->getPCU()->Add<long>(count[i]);

  if (m->getPCU()->Self()) return count;
  printf("%s Min: %f, Max: %f\n", name.c_str(), inputMin, inputMax);
  for (int i = 0; i < nbins; ++i) printf("%d\n", count[i]);
  return count;
}

HistogramStats printHistogramStats(Adapt* a)
{
  std::vector<double> lengths;
  std::vector<double> qualities;
  ma::stats(a->mesh, a->input->sizeField, lengths, qualities, true);
  std::vector<int> qualityHist = printHistogramData("\nQualities:", qualities, 0, 1, a->mesh);
  std::vector<int> lengthHist = printHistogramData("\nLengths:", lengths, 0, MAXLENGTH+1, a->mesh);
  return HistogramStats(qualityHist, lengthHist);
}

//Compare two histograms using L2 (Euclidean) distance.
double histogramDistance(const std::vector<int>& hist1, const std::vector<int>& hist2, bool normalize) 
{
  if (hist1.size() != hist2.size()) throw std::invalid_argument("Histograms must be the same size");

  size_t n = hist1.size();
  std::vector<double> h1(hist1.begin(), hist1.end());
  std::vector<double> h2(hist2.begin(), hist2.end());

  if (normalize) {
    double sum1 = std::accumulate(h1.begin(), h1.end(), 0.0);
    double sum2 = std::accumulate(h2.begin(), h2.end(), 0.0);
    if (sum1 > 0) for (double& val : h1) val /= sum1;
    if (sum2 > 0) for (double& val : h2) val /= sum2;
  }

  double sum_sq_diff = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double diff = h1[i] - h2[i];
    sum_sq_diff += diff * diff;
  }
  return std::sqrt(sum_sq_diff);
}


}
