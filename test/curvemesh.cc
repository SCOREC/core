#include <ma.h>
#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <SimUtil.h>
#include <maCurveMesh.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>
#include <apfField.h>
#include <maAdapt.h>

static void testInterpolationError(ma::Mesh* m, int entityDim,
    apf::DynamicVector & errors){
  ma::Iterator* it = m->begin(entityDim);
  ma::Entity* e;
  int id = 0;
  while ((e = m->iterate(it))) {
    ma::Model* g = m->toModel(e);
    if (m->getModelType(g) == m->getDimension())
      errors[id] = -1.0;
    else
      errors[id] = ma::interpolationError(m,e,11);
    id++;
  }
  m->end(it);
}

static void testElementSize(ma::Mesh* m)
{
  int dim = m->getDimension();
  double sizes[3] = {0,0,0};
  for( int d = 1; d <= dim; ++d){
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
      sizes[d-1] += v;
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
  printf("Total sizes for order %d %f %f %f\n",
    m->getCoordinateField()->getShape()->getOrder(),
    sizes[0],sizes[1],sizes[2]);
}

static void testBezier(const char* modelFile, const char* meshFile,
    const char* outFile, const int ne, const int nf)
{

  apf::DynamicMatrix edgeErrors(ne,6);
  apf::DynamicMatrix faceErrors(nf,6);

  for(int order = 1; order <= 6; ++order){
    apf::Mesh2* m2 = apf::loadMdsMesh(modelFile,meshFile);
    ma::Input* in = ma::configureIdentity(m2);
    ma::Adapt* adapt = new ma::Adapt(in);
    ma::BezierCurver bc(adapt,order);
    bc.run();
    ma::writePointSet(m2,1,21,outFile);
    ma::writePointSet(m2,2,21,outFile);

    testElementSize(m2);
    apf::DynamicVector ee(ne);
    apf::DynamicVector fe(nf);
    testInterpolationError(m2,1,ee);
    testInterpolationError(m2,2,fe);
    edgeErrors.setColumn(order-1,ee);
    faceErrors.setColumn(order-1,fe);

    m2->destroyNative();
    apf::destroyMesh(m2);
  }

  for(int order = 1; order <= 5; ++order)
    for(int id = 0; id < ne; ++id){
      if(edgeErrors(id,0) > 1.e-10)
        assert(edgeErrors(id,order-1) - edgeErrors(id,order) > 0.);
    }

  for(int order = 1; order <= 5; ++order)
    for(int id = 0; id < nf; ++id){
      if(faceErrors(id,0) > 1.e-10)
        assert(faceErrors(id,order-1) - faceErrors(id,order) > 0.);
    }
}

static void testGregory(const char* modelFile, const char* meshFile,
    const char* outFile, const int ne, const int nf)
{

    apf::Mesh2* m2 = apf::loadMdsMesh(modelFile,meshFile);
    ma::Input* in = ma::configureIdentity(m2);
    ma::Adapt* adapt = new ma::Adapt(in);
    ma::GregoryCurver gc(adapt,4);
    gc.run();
    testElementSize(m2);

    apf::DynamicVector ee(ne);
    apf::DynamicVector fe(nf);
    testInterpolationError(m2,1,ee);
    testInterpolationError(m2,2,fe);
    ma::writePointSet(m2,1,21,outFile);
    ma::writePointSet(m2,2,21,outFile);

    m2->destroyNative();
    apf::destroyMesh(m2);
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

  testBezier(modelFile,meshFile,outFile,ne,nf);
  testGregory(modelFile,meshFile,outFile,ne,nf);

  PCU_Comm_Free();
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  MPI_Finalize();
}
