#include <ma.h>
#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <SimUtil.h>
#include <maCurveMesh.h>
#include <apfField.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>

void testInterpolatedPoints(ma::Mesh* m){

  int dim = m->getDimension();
  for( int d = 1; d < dim; ++d){
    ma::Iterator* it = m->begin(d);
    ma::Entity* e;
    while ((e = m->iterate(it))) {
      ma::Vector samplePt, maxPt;
      double error = ma::interpolationErrorAtNodeXi(m,e,0,samplePt,maxPt);
      if (error > 1.e-7) {
        std::cerr << apf::Mesh::typeName[m->getType(e)] <<
            " is not being interpolated correctly by " <<
            m->getCoordinateField()->getShape()->getName() << "\n";
        abort();
      }
    }
    m->end(it);
  }
}
void testInterpolationError(ma::Mesh* m, int entityDim,
    apf::DynamicVector & errors){
  ma::Iterator* it = m->begin(entityDim);
  ma::Entity* e;
  int id = 0;
  while ((e = m->iterate(it))) {
    errors[id] = 0.0;
    ma::Model* g = m->toModel(e);
    if (m->getModelType(g) == m->getDimension()) continue;
    ma::Vector samplePt, maxPt;
    errors[id++]= ma::interpolationError(m,e,21,samplePt,maxPt);
  }
  m->end(it);
}
void testElementSize(ma::Mesh* m)
{
  int dim = m->getDimension();
  for( int d = 1; d < dim; ++d){
    ma::Iterator* it = m->begin(d);
    ma::Entity* e;
    while ((e = m->iterate(it))) {
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      if(v < 0){
        std::stringstream ss;
        ss << "error: " << apf::Mesh::typeName[m->getType(e)]
           << " size " << v
           << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
        abort();
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
}
int main(int argc, char** argv)
{
  assert(argc==4);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* outFile = argv[3];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  int ne = m->count(1);
  int nf = m->count(2);
  m->destroyNative();
  apf::destroyMesh(m);

  apf::DynamicMatrix edgeErrors(ne,6);
  apf::DynamicMatrix faceErrors(nf,6);

  for(int order = 1; order < 7; ++order){
    apf::Mesh2* m2 = apf::loadMdsMesh(modelFile,meshFile);
    ma::curveMeshToBezier(m2,order);
    testInterpolatedPoints(m2);
    testElementSize(m2);
    apf::DynamicVector ee(ne);
    apf::DynamicVector fe(nf);
    testInterpolationError(m2,1,ee);
    testInterpolationError(m2,2,fe);
    edgeErrors.setColumn(order-1,ee);
    faceErrors.setColumn(order-1,fe);
    ma::writePointSet(m2,1,21,outFile);
    ma::writePointSet(m2,2,21,outFile);
    m2->destroyNative();
    apf::destroyMesh(m2);
  }

  for(int id = 0; id < ne; ++id){
    apf::DynamicVector ee(6);
    edgeErrors.getRow(id,ee);
    if(ee[0] > 1e-12)
      printf("edge %d %.4e %.4e %.4e %.4e %.4e %.4e \n",
          id,ee[0],ee[1],ee[2],ee[3],ee[4],ee[5]);
  }
  for(int id = 0; id < nf; ++id){
    apf::DynamicVector fe(6);
    faceErrors.getRow(id,fe);
    if(fe[0] > 1e-12)
      printf("face %d %.4e %.4e %.4e %.4e %.4e %.4e \n",
          id,fe[0],fe[1],fe[2],fe[3],fe[4],fe[5]);
  }

  PCU_Comm_Free();
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  MPI_Finalize();
}

