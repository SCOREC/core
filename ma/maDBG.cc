/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maDBG.h"
#include "maShape.h"
#include <gmi.h>
#include <gmi_null.h>
#include <PCU.h>
#include <apf.h>
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

void addTargetLocation(ma::Adapt* a,
    const char* fieldName)
{
  ma::Mesh* m = a->mesh;
  apf::Field* paramField;
  paramField = m->findField(fieldName);
  if (paramField)
    apf::destroyField(paramField);

  paramField = apf::createFieldOn(m, fieldName, apf::VECTOR);
  ma::Entity* ent;
  ma::Iterator* it;
  it = m->begin(0);
  while ( (ent = m->iterate(it)) ){
    ma::Vector p;
    m->getParam(ent, p);
    ma::Vector x;
    if (m->getModelType(m->toModel(ent)) != 3)
      m->snapToModel(m->toModel(ent), p, x);
    else
      x = ma::getPosition(m, ent);

    x = x - ma::getPosition(m, ent);
    apf::setVector(paramField , ent, 0, x);
  }
  m->end(it);
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
  ma::getLinearQualitiesInMetricSpace(a->mesh, a->sizeField, lq_metric);
  ma::getLinearQualitiesInPhysicalSpace(a->mesh, lq_no_metric);

  colorEntitiesOfDimWithValues(a, a->mesh->getDimension(), lq_metric, "qual_metric");
  colorEntitiesOfDimWithValues(a, a->mesh->getDimension(), lq_no_metric, "qual_no_metric");

  // for snap debug
  if (a->mesh->canSnap())
    addTargetLocation(a, "target_for_snap");


  // setup file name and write the mesh
  std::stringstream ss;

  ss << a->input->debugFolder << "/";
  ss << std::setfill('0') << std::setw(3) << iter << "_";
  ss << prefix;

  writeMesh(a->mesh, ss.str().c_str(), "");

  apf::Field* colorField;
  colorField = a->mesh->findField("qual_metric");
  if (colorField)
    apf::destroyField(colorField);

  colorField = a->mesh->findField("qual_no_metric");
  if (colorField)
    apf::destroyField(colorField);

  colorField = a->mesh->findField("target_for_snap");
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

void createCavityMesh(ma::Adapt* a,
    ma::EntityArray& tets,
    const char* prefix)
{
  ma::Mesh* m = a->mesh;

  gmi_register_null();
  ma::Mesh* cavityMesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);

  size_t n = tets.getSize();
  for (size_t i = 0; i < n; ++i) {
    ma::Entity* down[4];
    ma::Entity* newVerts[4];
    m->getDownward(tets[i], 0, down);
    for (int j = 0; j < 4; j++) {
      ma::Vector position;
      m->getPoint(down[j], 0, position);
      newVerts[j] = cavityMesh->createVertex(NULL, position, ma::Vector(0,0,0));
    }
    apf::buildElement(cavityMesh, NULL, apf::Mesh::TET, newVerts);
  }

  cavityMesh->acceptChanges();
  std::stringstream ss;
  ss << a->input->debugFolder << "/";
  ss << prefix;


  apf::writeVtkFiles(ss.str().c_str(),cavityMesh);
  cavityMesh->destroyNative();
  apf::destroyMesh(cavityMesh);
}

void createCavityMesh(ma::Adapt* a,
    ma::EntitySet& tets,
    const char* prefix)
{
  ma::EntityArray tetsArray;
  tetsArray.setSize(tets.size());
  int count = 0;
  APF_ITERATE(ma::EntitySet,tets,it)
    tetsArray[count++] = *it;
  createCavityMesh(a, tetsArray, prefix);
}


}
