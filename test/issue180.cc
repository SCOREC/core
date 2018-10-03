#include <ma.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <parma.h>
#include <cassert>

void checktag(apf::Mesh* m, apf::MeshTag* t, double refVal) {
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  while( (e = m->iterate(it)) ) {
    assert( m->hasTag(e,t) );
    double v = 0;
    m->getDoubleTag(e, t, &v);
    assert(v==refVal);
  }
  m->end(it);
}

apf::Mesh2* createMesh2D() {
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false);
  apf::MeshEntity* v[3], *edges[3];
  apf::Vector3 points2D[3] =
  {apf::Vector3(0,0,0),apf::Vector3(1,0,0),apf::Vector3(1,1,0)};

  for (int i = 0; i < 3; ++i){
    v[i] = m->createVertex(m->findModelEntity(0,i),points2D[i],points2D[i]);
  }
  for (int i = 0; i < 3; ++i){
    apf::ModelEntity* edge = m->findModelEntity(1,i);
    apf::MeshEntity* ved[2] = {v[i],v[(i+1) % 3]};
    edges[i] = m->createEntity(apf::Mesh::EDGE,edge,ved);
  }
  apf::ModelEntity* faceModel = m->findModelEntity(2,0);
  m->createEntity(apf::Mesh::TRIANGLE,faceModel,edges);

  m->acceptChanges();
  m->verify();
  return m;
}

int main(int argc, char** argv) 
{ 
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  gmi_register_null();

  // create a mesh with a single triangle
  apf::Mesh2* m = createMesh2D();

  // create and set the tag
  apf::MeshTag* t = m->createDoubleTag("testtag2",1);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(2);
  double tagVal = 2;
  while( (e = m->iterate(it)) )
    m->setDoubleTag(e, t, &tagVal);
  m->end(it);

  // sanity check the tags
  checktag(m,t,tagVal);

  // uniform refinement
  ma::Input* in = ma::configureUniformRefine(m, 1);
  if (in->shouldSnap) {
    in->shouldSnap = false;
    PCU_ALWAYS_ASSERT(in->shouldTransferParametric);
  }
  in->shouldFixShape = false;
  ma::adapt(in);

  // sanity check the tags
  checktag(m,t,tagVal);

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
