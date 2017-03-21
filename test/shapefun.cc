#include <apfShape.h>
#include <apfMesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <pcu_util.h>

#include <iostream>

void testNodeValues(apf::EntityShape* shp, apf::Vector3 const* nodes, int nnodes)
{
  PCU_ALWAYS_ASSERT( shp->countNodes() == nnodes );
  double tol = 1.0e-15;
  for (int i = 0; i < nnodes; ++i) {
    apf::NewArray<double> values;
    shp->getValues(0, 0, nodes[i], values);
    for (int j = 0; j < nnodes; ++j) {
      if (j == i)
        PCU_ALWAYS_ASSERT(fabs(values[j] - 1) < tol);
      else
        PCU_ALWAYS_ASSERT(fabs(values[j] - 0) < tol);
    }
  }
}

void testP1LineNodeValues()
{
  apf::Vector3 nodes[2] = {
    apf::Vector3(-1,0,0),
    apf::Vector3(1,0,0) };
  apf::EntityShape* shp =
    apf::getLagrange(1)->getEntityShape(apf::Mesh::EDGE);
  testNodeValues(shp, nodes, 2);
}

void testP2LineNodeValues()
{
  apf::Vector3 nodes[3] = {
    apf::Vector3(-1,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,0,0) };
  apf::EntityShape* shp =
    apf::getLagrange(2)->getEntityShape(apf::Mesh::EDGE);
  testNodeValues(shp, nodes, 3);
}

void testP3LineNodeValues()
{
  apf::Vector3 nodes[4] = {
    apf::Vector3(-1,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(-1./3.,0,0),
    apf::Vector3(1./3.,0,0) };
  apf::EntityShape* shp =
    apf::getLagrange(3)->getEntityShape(apf::Mesh::EDGE);
  testNodeValues(shp, nodes, 4);
}

void testP1TriNodeValues() {
  apf::Vector3 nodes[3] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0) };
  apf::EntityShape* shp =
    apf::getLagrange(1)->getEntityShape(apf::Mesh::TRIANGLE);
  testNodeValues(shp, nodes, 3);
}

void testP2TriNodeValues() {
  apf::Vector3 nodes[6] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(.5,0,0),
    apf::Vector3(.5,.5,0),
    apf::Vector3(0,.5,0) };
  apf::EntityShape* shp =
    apf::getLagrange(2)->getEntityShape(apf::Mesh::TRIANGLE);
  testNodeValues(shp, nodes, 6);
}

void testP3TriNodeValues() {
  apf::Vector3 nodes[10] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(1./3.,0,0),
    apf::Vector3(2./3.,0,0),
    apf::Vector3(2./3.,1./3.,0),
    apf::Vector3(1./3.,2./3.,0),
    apf::Vector3(0,2./3.,0),
    apf::Vector3(0,1./3.,0),
    apf::Vector3(1./3.,1./3.,0) };
  apf::EntityShape* shp =
    apf::getLagrange(3)->getEntityShape(apf::Mesh::TRIANGLE);
  testNodeValues(shp, nodes, 10);
}

void testP1TetNodeValues() {
  apf::Vector3 nodes[4] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(0,0,1) };
    apf::EntityShape* shp =
      apf::getLagrange(1)->getEntityShape(apf::Mesh::TET);
    testNodeValues(shp, nodes, 4);
}

void testP2TetNodeValues() {
  apf::Vector3 nodes[10] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(0,0,1),
    apf::Vector3(0.5,0,0),
    apf::Vector3(0.5,0.5,0),
    apf::Vector3(0,0.5,0),
    apf::Vector3(0,0,0.5),
    apf::Vector3(0.5,0,0.5),
    apf::Vector3(0,0.5,0.5) };
    apf::EntityShape* shp =
      apf::getLagrange(2)->getEntityShape(apf::Mesh::TET);
    testNodeValues(shp, nodes, 10);
}

void testP3TetNodeValues() {
  apf::Vector3 nodes[20] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(0,0,1),
    apf::Vector3(1./3.,0,0),
    apf::Vector3(2./3.,0,0),
    apf::Vector3(2./3.,1./3.,0),
    apf::Vector3(1./3.,2./3.,0),
    apf::Vector3(0,2./3.,0),
    apf::Vector3(0,1./3.,0),
    apf::Vector3(0,0,1./3.),
    apf::Vector3(0,0,2./3.),
    apf::Vector3(2./3.,0,1./3.),
    apf::Vector3(1./3.,0,2./3.),
    apf::Vector3(0,2./3.,1./3.),
    apf::Vector3(0,1./3.,2./3.),
    apf::Vector3(1./3.,1./3.,0),
    apf::Vector3(1./3.,0,1./3.),
    apf::Vector3(1./3,1./3.,1./3.),
    apf::Vector3(0,1./3.,1./3.) };
  apf::EntityShape* shp =
    apf::getLagrange(3)->getEntityShape(apf::Mesh::TET);
  testNodeValues(shp, nodes, 20);
}


void testPrismNodeValues()
{
  apf::Vector3 nodes[6];
  int i = 0;
  for (int xi2 = -1; xi2 < 2; xi2 += 2)
  for (int xi1 = 0; xi1 < 2; ++xi1)
  for (int xi0 = 0; xi0 < 2; ++xi0) {
    if (xi0 == 1 && xi1 == 1)
      continue;
    nodes[i] = apf::Vector3(xi0,xi1,xi2);
    ++i;
  }
  apf::EntityShape* shp = apf::getLagrange(1)->getEntityShape(apf::Mesh::PRISM);
  testNodeValues(shp, nodes, 6);
}

void testPyramidNodeValues()
{
  apf::Vector3 nodes[5] =
  {apf::Vector3(-1,-1,-1),
   apf::Vector3( 1,-1,-1),
   apf::Vector3( 1, 1,-1),
   apf::Vector3(-1, 1,-1),
   apf::Vector3( 0, 0, 1)};
  apf::EntityShape* shp = apf::getLagrange(1)->getEntityShape(apf::Mesh::PYRAMID);
  testNodeValues(shp, nodes, 5);
}

void testQuadrilateralNodeValues() {
  apf::Vector3 nodes[4] =
  {apf::Vector3(-1,-1, 0),
   apf::Vector3( 1,-1, 0),
   apf::Vector3( 1, 1, 0),
   apf::Vector3(-1, 1, 0)};
  apf::EntityShape* shp = apf::getLagrange(1)->getEntityShape(apf::Mesh::QUAD);
  testNodeValues(shp, nodes, 4);
  /*second order shape functions*/
  apf::Vector3 nodes2[9] =
  {apf::Vector3(-1,-1, 0),
   apf::Vector3( 1,-1, 0),
   apf::Vector3( 1, 1, 0),
   apf::Vector3(-1, 1, 0),
   apf::Vector3( 0,-1, 0),
   apf::Vector3( 1, 0, 0),
   apf::Vector3( 0, 1, 0),
   apf::Vector3(-1, 0, 0),
   apf::Vector3( 0, 0, 0)};

  apf::EntityShape* shp2 = apf::getLagrange(2)->getEntityShape(apf::Mesh::QUAD);
  testNodeValues(shp2, nodes2, 9);
}

void testVolume(int type, apf::Vector3 const* points, double volume)
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false);
  apf::MeshEntity* e = apf::buildOneElement(
      m, m->findModelEntity(3,0), type, points);
  double v = apf::measure(m,e);
  PCU_ALWAYS_ASSERT( fabs(v - volume) < 1e-10 );
  m->destroyNative();
  apf::destroyMesh(m);
}

void testPrismVolume()
{
  apf::Vector3 points[6] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,1),
   apf::Vector3(1,0,1),
   apf::Vector3(0,1,1)};
  testVolume(apf::Mesh::PRISM, points, 0.5);
}

void testPyramidVolume()
{
  apf::Vector3 points[5] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(1,1,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,3)};
  testVolume(apf::Mesh::PYRAMID, points, 1);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  testP1LineNodeValues();
  testP2LineNodeValues();
  testP3LineNodeValues();
  testP1TriNodeValues();
  testP2TriNodeValues();
  testP3TriNodeValues();
  testP1TetNodeValues();
  testP2TetNodeValues();
  testP3TetNodeValues();
  testPrismNodeValues();
  testPyramidNodeValues();
  testQuadrilateralNodeValues();
  testPrismVolume();
  testPyramidVolume();
  PCU_Comm_Free();
  MPI_Finalize();
}
