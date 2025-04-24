/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maDBG.h"
#include "maShape.h"
#include "maAdapt.h"
#include "maRefine.h"
#include <gmi.h>
#include <gmi_null.h>
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

static double PI = 3.14159265359;

namespace ma_dbg {


void writeMesh(ma::Mesh* m,
    const char* prefix,
    const char* suffix,
    int dim=-1)
{
  std::stringstream ss;
  if (std::string(suffix) != "")
    ss << prefix << "/" << suffix;
  else
    ss << prefix;
  std::string tmp = ss.str();
  const char* fileName = tmp.c_str();
  apf::writeVtkFiles(fileName, m, dim);
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

void addParamCoords(ma::Adapt* a,
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
    apf::setVector(paramField , ent, 0, p);
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


  if (dim == 0)
    colorField = apf::createFieldOn(m, fieldName, apf::SCALAR);
  else
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

  // parametric coordinates
  addParamCoords(a, "param_coords");

  // setup file name and write the mesh
  std::stringstream ss;

  if (a->input->debugFolder) {
    ss << a->input->debugFolder << "/";
  }
  ss << std::setfill('0') << std::setw(3) << iter << "_";
  ss << prefix;

  writeMesh(a->mesh, ss.str().c_str(), "", -1);

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

  colorField = a->mesh->findField("param_coords");
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

  colorEntitiesOfDimWithValues(a, dim, ent_flags, flagName);

  // setup file name and write the mesh
  std::stringstream ss;
  if (a->input->debugFolder) {
    ss << a->input->debugFolder << "/";
  }
  ss << prefix << "_" << std::setfill('0') << std::setw(3) << iter;

  writeMesh(a->mesh, ss.str().c_str(), "", dim);

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
  ma::Mesh* cavityMesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false, m->getPCU());

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
  if (a->input->debugFolder) {
    ss << a->input->debugFolder << "/";
  }
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

static apf::Vector3 getPointOnEllipsoid(
    apf::Vector3 center,
    apf::Vector3 abc,
    apf::Matrix3x3 rotation,
    double scaleFactor,
    double u,
    double v)
{
  apf::Vector3 result;
  result[0] = abc[0] * cos(u) * cos(v);
  result[1] = abc[1] * cos(u) * sin(v);
  result[2] = abc[2] * sin(u);

  result = result * scaleFactor;

  result = rotation * result + center;
  return result;
}

static void makeEllipsoid(
    apf::Mesh2* msf,
    apf::Mesh2* mesh,
    apf::Field* sizes,
    apf::Field* frames,
    apf::MeshEntity* ent,
    int node,
    double scaleFactor,
    int sampleSize[2])
{
  // first get the coordinate at node location
  apf::Vector3 xi;
  apf::Vector3 center;
  apf::FieldShape* fs = apf::getShape(sizes);
  fs->getNodeXi(mesh->getType(ent), node, xi);
  apf::MeshElement* me = apf::createMeshElement(mesh, ent);
  apf::mapLocalToGlobal(me, xi, center);
  apf::destroyMeshElement(me);


  // second get the sizes and frames at node
  apf::Vector3 abc;
  apf::getVector(sizes, ent, node, abc);

  apf::Matrix3x3 rotation;
  apf::getMatrix(frames, ent, node, rotation);


  double U0 = 0.0;
  double U1 = 2 * PI;
  double V0 = -PI/2.;
  double V1 =  PI/2.;
  int n = sampleSize[0];
  int m = sampleSize[1];
  double dU = (U1 - U0) / (n-1);
  double dV = (V1 - V0) / (m-1);

  // make the array of vertex coordinates in the physical space
  std::vector<ma::Vector> ps;
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      double u = U0 + i * dU;
      double v = V0 + j * dV;
      apf::Vector3 pt = getPointOnEllipsoid(center, abc, rotation, scaleFactor, u, v);
      ps.push_back(pt);
    }
  }
  // make the vertexes and set the coordinates using the array
  std::vector<apf::MeshEntity*> vs;
  for (size_t i = 0; i < ps.size(); i++) {
    apf::MeshEntity* newVert = msf->createVert(0);
    msf->setPoint(newVert, 0, ps[i]);
    vs.push_back(newVert);
  }

  PCU_ALWAYS_ASSERT(vs.size() == ps.size());

  apf::MeshEntity* v[3];
  // make the lower/upper t elems
  for (int i = 0; i < n-1; i++) {
    for (int j = 0; j < m-1; j++) {
      // upper triangle
      v[0] = vs[(i + 0) + n * (j + 0)];
      v[1] = vs[(i + 0) + n * (j + 1)];
      v[2] = vs[(i + 1) + n * (j + 0)];
      apf::buildElement(msf, 0, apf::Mesh::TRIANGLE, v);
      // upper triangle
      v[0] = vs[(i + 0) + n * (j + 1)];
      v[1] = vs[(i + 1) + n * (j + 1)];
      v[2] = vs[(i + 1) + n * (j + 0)];
      apf::buildElement(msf, 0, apf::Mesh::TRIANGLE, v);
    }
  }
}

void visualizeSizeField(
    apf::Mesh2* m,
    apf::Field* sizes,
    apf::Field* frames,
    int sampleSize[2],
    double userScale,
    const char* outputPrefix)
{
  // create the size-field visualization mesh
  apf::Mesh2* msf = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false, m->getPCU());

  apf::FieldShape* fs = apf::getShape(sizes);
  int dim = m->getDimension();

  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!fs->hasNodesIn(d)) continue;
    it = m->begin(d);
    while ( (ent = m->iterate(it)) ) {
      int type = m->getType(ent);
      int non = fs->countNodesOn(type);
      for (int n = 0; n < non; n++)
	makeEllipsoid(msf, m, sizes, frames, ent, n, userScale , sampleSize);
    }
    m->end(it);
  }

  apf::writeVtkFiles(outputPrefix, msf);

  msf->destroyNative();
  apf::destroyMesh(msf);
}

struct SplitByTag : public ma::Predicate
{
  SplitByTag(ma::Adapt* a_, int type_, int tag_) :
    a(a_), type(type_), tag(tag_) {}
  bool operator() (ma::Entity* e) {
    ma::Mesh* m = a->mesh;
    int mtype = m->getModelType(m->toModel(e));
    int mtag  = m->getModelTag(m->toModel(e));
    if (mtype != type) return false;
    if (mtag  != tag ) return false;
    return true;
  }
  ma::Adapt* a;
  int type;
  int tag;
};

void uniformAdaptByModelTag(
    apf::Mesh2* m,
    int mtype,
    int mtag,
    int level)
{
  ma::Input* in = ma::makeAdvanced(ma::configureUniformRefine(m, 0));
  ma::validateInput(in);
  ma::Adapt* a = new ma::Adapt(in);
  for(int i = 0; i < level; i++) {
    SplitByTag p(a, mtype, mtag);
    ma::markEntities(a, 1, p, ma::SPLIT, ma::NEED_NOT_SPLIT,
    	ma::DONT_SPLIT | ma::NEED_NOT_SPLIT);
    PCU_ALWAYS_ASSERT(ma::checkFlagConsistency(a,1,ma::SPLIT));
    ma::Refine* r = a->refine;
    ma::resetCollection(r);
    ma::collectForTransfer(r);
    ma::collectForMatching(r);
    ma::addAllMarkedEdges(r);
    ma::splitElements(r);
    ma::processNewElements(r);
    ma::destroySplitElements(r);
    ma::forgetNewEntities(r);
  }
  delete(a);
}

}
