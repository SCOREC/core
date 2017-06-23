#include <crv.h>

#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <pcu_util.h>

class Linear : public ma::IsotropicFunction
{
  public:
    Linear(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m,lower,upper);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double x = (p[0] - lower[0])/(upper[0] - lower[0]);
      return average*(4*x+2)/3;
    }
  private:
    ma::Mesh* mesh;
    double average;
    ma::Vector lower;
    ma::Vector upper;
};

static void testInterpolationError(apf::Mesh* m, int entityDim,
    apf::DynamicVector & errors){
  apf::MeshIterator* it = m->begin(entityDim);
  apf::MeshEntity* e;
  int id = 0;
  while ((e = m->iterate(it))) {
    apf::ModelEntity* g = m->toModel(e);
    if (m->getModelType(g) == m->getDimension())
      errors[id] = -1.0;
    else
      errors[id] = crv::interpolationError(m,e,11);
    id++;
  }
  m->end(it);
}

static void testElementSize(apf::Mesh* m)
{
  int dim = m->getDimension();
  double sizes[3] = {0,0,0};
  for( int d = 1; d <= dim; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      if(!m->isOwned(e)) continue;
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      if(v < 0){
        std::stringstream ss;
        ss << "error: " << apf::Mesh::typeName[m->getType(e)]
           << " size " << v
           << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
        PCU_ALWAYS_ASSERT(v < 0);
      }
      sizes[d-1] += v;
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }

  PCU_Add_Doubles(&(sizes[0]),3);
  if(!PCU_Comm_Self())
    printf("Total sizes for order %d %f %f %f\n",
      m->getShape()->getOrder(),
      sizes[0],sizes[1],sizes[2]);
}

static void testInterpolating(const char* modelFile, const char* meshFile,
    const int ne, const int nf)
{
  for(int order = 1; order <= 2; ++order){
    apf::Mesh2* m2 = apf::loadMdsMesh(modelFile,meshFile);
    apf::changeMeshShape(m2,apf::getLagrange(order),true);
    crv::InterpolatingCurver ic(m2,order);
    ic.run();
    testElementSize(m2);
    apf::DynamicVector ee(ne);
    apf::DynamicVector fe(nf);
    testInterpolationError(m2,1,ee);
    testInterpolationError(m2,2,fe);
    m2->destroyNative();
    apf::destroyMesh(m2);
  }
}

static void testBezier(const char* modelFile, const char* meshFile,
    const int ne, const int nf)
{

  apf::DynamicMatrix edgeErrors(ne,6);
  apf::DynamicMatrix faceErrors(nf,6);
  // check blended shapes + interpolation error
  // interpolation error decreases as order is increased
  for(int order = 1; order <= 6; ++order){
    apf::Mesh2* m2 = apf::loadMdsMesh(modelFile,meshFile);
    crv::BezierCurver bc(m2,order,2);
    bc.run();
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
        PCU_ALWAYS_ASSERT(edgeErrors(id,order-1) - edgeErrors(id,order) > 0.);
    }

  for(int order = 1; order <= 5; ++order)
    for(int id = 0; id < nf; ++id){
      if(faceErrors(id,0) > 1.e-10)
        PCU_ALWAYS_ASSERT(faceErrors(id,order-1) - faceErrors(id,order) > 0.);
    }
  // check some refinement
  for(int order = 2; order <= 4; ++order){
    apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
    crv::BezierCurver bc(m,order,0);
    bc.run();
    Linear sf(m);
    ma::Input* in = ma::configure(m,&sf);
    in->shouldSnap = true;
    in->shouldTransferParametric = true;
    in->maximumIterations = 1;
    crv::adapt(in);
    m->destroyNative();
    apf::destroyMesh(m);
  }
  // check some refinement for blended shapes
  for(int order = 2; order <= 4; ++order){
    apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
    crv::BezierCurver bc(m,order,1);
    bc.run();
    Linear sf(m);
    ma::Input* in = ma::configure(m,&sf);
    in->shouldSnap = true;
    in->shouldTransferParametric = true;
    in->maximumIterations = 1;
    crv::adapt(in);
    m->destroyNative();
    apf::destroyMesh(m);
  }

}

static void testGregory(const char* modelFile, const char* meshFile,
    const int ne, const int nf)
{
  for(int order = 0; order <= 2; ++order){
    apf::Mesh2* m2 = apf::loadMdsMesh(modelFile,meshFile);
    crv::GregoryCurver gc(m2,4,order);
    gc.run();
    testElementSize(m2);
    apf::DynamicVector ee(ne);
    apf::DynamicVector fe(nf);
    testInterpolationError(m2,1,ee);
    testInterpolationError(m2,2,fe);
    m2->destroyNative();
    apf::destroyMesh(m2);
  }
}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  int ne = m->count(1);
  int nf = m->count(2);
  m->destroyNative();
  apf::destroyMesh(m);

  testInterpolating(modelFile,meshFile,ne,nf);
  testBezier(modelFile,meshFile,ne,nf);
  testGregory(modelFile,meshFile,ne,nf);

  PCU_Comm_Free();
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
  MPI_Finalize();
}
