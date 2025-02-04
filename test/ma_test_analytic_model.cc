#include <crv.h>
#include <crvAdapt.h>
#include <crvBezier.h>
#include <crvBezierShapes.h>
#include <crvQuality.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <lionPrint.h>

#include <math.h>
#include <pcu_util.h>

const double pi = 3.14159265359;

class myCallback : public apf::BuildCallback
{
  public:
    myCallback(apf::Mesh2* inMesh) : mesh(inMesh) {};
    ~myCallback() {};
    virtual void call(apf::MeshEntity* e)
    {
      int eType = mesh->getType(e);
      int mType = mesh->getModelType(mesh->toModel(e));
      int mTag  = mesh->getModelTag(mesh->toModel(e));
      printf("entity type is %d\n", eType);
      printf("model type/tag is %d/%d\n", mType, mTag);
    }
  private:
    apf::Mesh2* mesh;
};

void face0(double const p[2], double x[3], void*)
{
  x[0] = cos(p[0]) * cos(p[1]);
  x[1] = sin(p[0]) * cos(p[1]);
  x[2] = sin(p[1]);
}
agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

gmi_model* makeSphere()
{
  gmi_model* model = gmi_make_analytic();
  int faPer[2] = {1, 0};
  double faRan[2][2] = {{0,6.28318530718},{-1.57079632679,1.57079632679}};
  /* gmi_ent* f = gmi_add_analytic(model, 2, 0, face0, faPer, faRan, 0); */
  gmi_add_analytic(model, 2, 0, face0, faPer, faRan, 0);

  gmi_add_analytic_region(model, 1);
  return model;
}

apf::Mesh2* createSphereMesh(pcu::PCU *PCUObj)
{
  gmi_model* model = makeSphere();
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false, PCUObj);
  apf::MeshEntity* allVs[5];
  apf::Vector3 p0( cos(0.), sin(0.),  0.);
  apf::Vector3 p1( cos(2.*pi/3.), sin(2.*pi/3.),  0.);
  apf::Vector3 p2( cos(4.*pi/3.), sin(4.*pi/3.),  0.);
  apf::Vector3 p3( 0.,  0.,  1.);
  apf::Vector3 p4( 0.,  0., -1.);


  apf::Vector3 uv0( 0.,  0.,  0.);
  apf::Vector3 uv1( 2.*pi/3.,  0.,  0.);
  apf::Vector3 uv2( 4.*pi/3.,  0.,  0.);
  apf::Vector3 uv3( 0.,  pi/2.,  0.);
  apf::Vector3 uv4( 0., -pi/2.,  0.);

  allVs[0] = m->createVertex(m->findModelEntity(2,0),p0,uv0);
  allVs[1] = m->createVertex(m->findModelEntity(2,0),p1,uv1);
  allVs[2] = m->createVertex(m->findModelEntity(2,0),p2,uv2);
  allVs[3] = m->createVertex(m->findModelEntity(2,0),p3,uv3);
  allVs[4] = m->createVertex(m->findModelEntity(2,0),p4,uv4);


  // create edges
  apf::MeshEntity* evs[2];
  apf::MeshEntity* es[3];

  evs[0] = allVs[0];
  evs[1] = allVs[1];
  es[0] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::EDGE, evs);

  evs[0] = allVs[1];
  evs[1] = allVs[2];
  es[1] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::EDGE, evs);

  evs[0] = allVs[2];
  evs[1] = allVs[0];
  es[2] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::EDGE, evs);


  apf::MeshEntity* sharedFace = m->createEntity(apf::Mesh::TRIANGLE, m->findModelEntity(3,1), es);


  apf::MeshEntity* fvs[3];
  apf::MeshEntity* fs[4];


  apf::MeshEntity* vs[4] = {allVs[0],allVs[1],allVs[2],allVs[3]};
  fs[0] = sharedFace;

  fvs[0] = vs[apf::tet_tri_verts[1][0]];
  fvs[1] = vs[apf::tet_tri_verts[1][1]];
  fvs[2] = vs[apf::tet_tri_verts[1][2]];
  fs[1] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, fvs);

  fvs[0] = vs[apf::tet_tri_verts[2][0]];
  fvs[1] = vs[apf::tet_tri_verts[2][1]];
  fvs[2] = vs[apf::tet_tri_verts[2][2]];
  fs[2] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, fvs);


  fvs[0] = vs[apf::tet_tri_verts[3][0]];
  fvs[1] = vs[apf::tet_tri_verts[3][1]];
  fvs[2] = vs[apf::tet_tri_verts[3][2]];
  fs[3] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, fvs);

  m->createEntity(apf::Mesh::TET, m->findModelEntity(3,1), fs);


  vs[0] = allVs[0];
  vs[1] = allVs[1];
  vs[2] = allVs[4];
  vs[3] = allVs[2];

  fvs[0] = vs[apf::tet_tri_verts[0][0]];
  fvs[1] = vs[apf::tet_tri_verts[0][1]];
  fvs[2] = vs[apf::tet_tri_verts[0][2]];
  fs[0] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, fvs);

  fs[1] = sharedFace;

  fvs[0] = vs[apf::tet_tri_verts[2][0]];
  fvs[1] = vs[apf::tet_tri_verts[2][1]];
  fvs[2] = vs[apf::tet_tri_verts[2][2]];
  fs[2] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, fvs);


  fvs[0] = vs[apf::tet_tri_verts[3][0]];
  fvs[1] = vs[apf::tet_tri_verts[3][1]];
  fvs[2] = vs[apf::tet_tri_verts[3][2]];
  fs[3] = apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, fvs);


  m->createEntity(apf::Mesh::TET, m->findModelEntity(3,1), fs);

  m->acceptChanges();
  m->verify();
  return m;
}


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);

  apf::Mesh2* m = createSphereMesh(&PCUObj);
  m->verify();

  apf::writeVtkFiles("initial_mesh_on_analytic_model", m);

  const ma::Input* in = ma::configureUniformRefine(m, 2);
  ma::adapt(in);

  apf::writeVtkFiles("adapted_mesh_on_analytic_model", m);

  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}
