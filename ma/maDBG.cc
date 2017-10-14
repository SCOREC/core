/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maDBG.h"
#include "maShape.h"
#include <PCU.h>
#include <apf.h>
/* #include <apfMesh2.h> */
#include <apfMDS.h>
#include <apfConvert.h>
#include <apfDynamicArray.h>
#include <pcu_util.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <cassert>
#include <stdarg.h>


namespace ma_dbg {


void writeMesh(ma::Mesh* m,
    const char* prefix,
    const char* suffix)
{
  std::stringstream ss;
  if (std::string(suffix) != "")
    ss << prefix << "/" << suffix;
  else
    ss << prefix;
  std::string tmp = ss.str();
  const char* fileName = tmp.c_str();
  apf::writeVtkFiles(fileName, m);
}

void colorEntitiesOfDimWithValues(ma::Adapt* a,
    int dim,
    const std::vector<double> & vals,
    const char* fieldName)
{
  ma::Mesh* m = a->mesh;
  apf::Field* colorField;
  colorField = m->findField(fieldName);
  if (colorField)
    apf::destroyField(colorField);


  colorField = apf::createField(m, fieldName, apf::SCALAR, apf::getConstant(dim));
  ma::Entity* ent;
  ma::Iterator* it;
  it = m->begin(dim);
  int index = 0;
  while ( (ent = m->iterate(it)) ){
    double color = (double) vals[index];
    apf::setComponents(colorField , ent, 0, &color);
    index++;
  }
  m->end(it);
}

void measureLinearQualties(ma::Adapt* a,
    std::vector<double> &lq,
    bool inMetric)
{
  ma::Mesh* m = a->mesh;
  int dim = m->getDimension();
  ma::SizeField* sf = a->sizeField;
  ma::Entity* e;
  ma::Iterator* it;

  it = m->begin(dim);
  while ( (e = m->iterate(it)) ) {
    if (inMetric) {
      lq.push_back(ma::measureElementQuality(m, sf, e, true));
    }
    else {
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
    double q;
    if (A < 0)
      q = -48. * (A*A) / (s*s);
    else
      q =  48. * (A*A) / (s*s);
    lq.push_back(q);
  }
  else if (m->getType(e) == apf::Mesh::TET){
    ma::Vector p[4];
    ma::getVertPoints(m, e, p);
    lq.push_back(ma::measureLinearTetQuality(p));
  }
    }
  }
  m->end(it);
}

void evaluateFlags(ma::Adapt* a,
    int dim,
    int flag,
    std::vector<double> &flgs)
{
  ma::Mesh* m = a->mesh;
  ma::Entity* e;
  ma::Iterator* it;

  it = m->begin(dim);
  while ( (e = m->iterate(it)) ) {
    bool hasFlag = ma::getFlag(a, e, flag);
    flgs.push_back(hasFlag ? 1.0 : 0.0);
  }
  m->end(it);
}

void dumpMeshWithQualities(ma::Adapt* a,
    int iter,
    const char* prefix)
{
  // measure qualities
  std::vector<double> lq_metric;
  std::vector<double> lq_no_metric;
  measureLinearQualties(a, lq_metric, true);
  measureLinearQualties(a, lq_no_metric, false);


  for (size_t i = 0; i < lq_metric.size(); i++) {
    lq_metric[i] = cbrt(lq_metric[i]);
    lq_no_metric[i] = cbrt(lq_no_metric[i]);
  }
  colorEntitiesOfDimWithValues(a, a->mesh->getDimension(), lq_metric, "qual_metric");
  colorEntitiesOfDimWithValues(a, a->mesh->getDimension(), lq_no_metric, "qual_no_metric");

  // setup file name and write the mesh
  std::stringstream ss;
  ss << a->input->debugFolder << "/";
  ss << prefix << "_" << std::setfill('0') << std::setw(3) << iter;

  writeMesh(a->mesh, ss.str().c_str(), "");

  apf::Field* colorField;
  colorField = a->mesh->findField("qual_metric");
  if (colorField)
    apf::destroyField(colorField);

  colorField = a->mesh->findField("qual_no_metric");
  if (colorField)
    apf::destroyField(colorField);
}

void dumpMeshWithFlag(ma::Adapt* a,
    int iter,
    int dim,
    int flag,
    const char* flagName,
    const char* prefix)
{
  // evaluate flags
  std::vector<double> ent_flags;
  evaluateFlags(a, dim, flag, ent_flags);

  colorEntitiesOfDimWithValues(a, a->mesh->getDimension(), ent_flags, flagName);

  // setup file name and write the mesh
  std::stringstream ss;
  ss << a->input->debugFolder << "/";
  ss << prefix << "_" << std::setfill('0') << std::setw(3) << iter;

  writeMesh(a->mesh, ss.str().c_str(), "");

  apf::Field* colorField;
  colorField = a->mesh->findField(flagName);
  if (colorField)
    apf::destroyField(colorField);
}


}
