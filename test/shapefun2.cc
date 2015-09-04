#include <apfShape.h>
#include <apfMesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <crv.h>
#include <PCU.h>
#include <cassert>

/* Test all shape functions by comparing them to the base linear,
 * by changing the mesh shape to the new shape functions
 * This assumes the base linear is correct, which is not tested here
 *
 * Only holds for interpolating shapes,
 * where the default projection is correct
 *
 */

// This a collection of coordinates for each shape
static apf::Vector3 const edge[2] = {
  apf::Vector3(-1,0,0),
  apf::Vector3(1,0,0)
};
static apf::Vector3 const tri[3] = {
  apf::Vector3(0,0,0),
  apf::Vector3(1,0,0),
  apf::Vector3(0,1,0)
};
static apf::Vector3 const quad[4] = {
  apf::Vector3(-1,-1,0),
  apf::Vector3(1,-1,0),
  apf::Vector3(1,1,0),
  apf::Vector3(-1,1,0)
};
static apf::Vector3 const tet[4] = {
  apf::Vector3(0,0,0),
  apf::Vector3(1,0,0),
  apf::Vector3(0,1,0),
  apf::Vector3(0,0,1)
};
static apf::Vector3 const hex[8] = {
  apf::Vector3(-1,-1,-1),
  apf::Vector3(1,-1,-1),
  apf::Vector3(1,1,-1),
  apf::Vector3(-1,1,-1),
  apf::Vector3(-1,-1,1),
  apf::Vector3(1,-1,1),
  apf::Vector3(1,1,1),
  apf::Vector3(-1,1,1)
};
static apf::Vector3 const prism[6] = {
  apf::Vector3(0,0,-1),
  apf::Vector3(1, 0, -1),
  apf::Vector3(0, 1, -1),
  apf::Vector3(0, 0, 1),
  apf::Vector3(1, 0, 1),
  apf::Vector3(0, 1, 1)
};
static apf::Vector3 const pyramid[5] = {
  apf::Vector3(-1,-1,-1),
  apf::Vector3( 1,-1,-1),
  apf::Vector3( 1, 1,-1),
  apf::Vector3(-1, 1,-1),
  apf::Vector3( 0, 0, 1)
};

static apf::Vector3 const* const points[apf::Mesh::TYPES] =
{0,edge,tri,quad,tet,hex,prism,pyramid};


/* Take a linear mesh and a higher-order mesh, and make sure
 * on a reference element, they give identical results
 *
 * This ensures that on a linear element, the shape functions being
 * tested reduce down exactly to the linear mesh
 *
 */
static void checkEntityShape(apf::Mesh* lm, apf::Mesh* m,
    apf::MeshEntity* le, apf::MeshEntity* e)
{


  apf::Vector3 lvalue, value;
  apf::Matrix3x3 lJacobian, Jacobian;

  int type = lm->getType(le);
  assert(lm->getType(le) == m->getType(e));

  int typeDim = apf::Mesh::typeDimension[type];

  apf::MeshElement* lme = apf::createMeshElement(lm,le);
  apf::MeshElement* me = apf::createMeshElement(m,e);

  apf::Element* lelem = apf::createElement(lm->getCoordinateField(),le);
  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);

  // xi's in each dimension to test things at,
  // chosen to avoid hitting a nodeXi exactly,
  // two points per dimension are chosen, one on the interior,
  // one on the boundary
  static apf::Vector3 xi[2][3] = {{apf::Vector3(1./3.,0,0),
      apf::Vector3(1./4,1./3,0), apf::Vector3(1./6,1./4,1./2)},
      {apf::Vector3(1.,0,0), apf::Vector3(1./6,0,0),
          apf::Vector3(0,1./4,1./3)}};
  for (int i = 0; i < 2; ++i){
    apf::getVector(lelem,xi[i][typeDim-1],lvalue);
    apf::getVector(elem,xi[i][typeDim-1],value);

    apf::getJacobian(lme,xi[i][typeDim-1],lJacobian);
    apf::getJacobian(me,xi[i][typeDim-1],Jacobian);
    // check particular point location
    assert(std::abs((lvalue-value).getLength()) < 1e-15);
    // check length/area/volume
    assert(std::abs(apf::measure(lme)-apf::measure(me)) < 1e-15);

    // check Jacobians are identical
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        assert(std::abs(Jacobian[i][j]-lJacobian[i][j]) < 1e-15);
  }
  // clean up
  apf::destroyMeshElement(lme);
  apf::destroyMeshElement(me);

  apf::destroyElement(lelem);
  apf::destroyElement(elem);
}
/* For every shape function we test, for every type we test,
 * create two meshes
 */
static void checkFieldShape(apf::FieldShape* fs)
{
  for (int type = 1; type < apf::Mesh::TYPES; ++type){
    apf::EntityShape* es = fs->getEntityShape(type);
    if(es){
      int typeDim = apf::Mesh::typeDimension[type];

      gmi_model* lmodel = gmi_load(".null");
      gmi_model* model = gmi_load(".null");

      apf::Mesh2* lm = apf::makeEmptyMdsMesh(lmodel, typeDim, false);
      apf::Mesh2* m = apf::makeEmptyMdsMesh(model, typeDim, false);

      apf::MeshEntity* le = apf::buildOneElement(
          lm, lm->findModelEntity(typeDim,0), type, points[type]);
      apf::MeshEntity* e = apf::buildOneElement(
          m, m->findModelEntity(typeDim,0), type, points[type]);

      lm->acceptChanges();
      m->acceptChanges();

      // change one mesh to the new shape functions
      apf::changeMeshShape(m,fs,true);

      checkEntityShape(lm,m,le,e);

      lm->destroyNative();
      m->destroyNative();

      apf::destroyMesh(lm);
      apf::destroyMesh(m);
    }
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();

  // put fieldShapes to test here
  apf::FieldShape* fs[4] = {apf::getLagrange(2),apf::getSerendipity(),
  crv::getBezier(3,1),crv::getBezier(3,2)};

  for (int i = 0; i < 4; ++i)
    checkFieldShape(fs[i]);

  PCU_Comm_Free();
  MPI_Finalize();
}

