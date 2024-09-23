#include <iostream>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <mpi.h>
#include <lionPrint.h>
#include <apf.h>
#include <apfShape.h>
#include <pcu_util.h>
#include <gmi_mesh.h>

/*
 * This test sets a nodal point matrix on a linear mesh. It then computes the gradient 
 *  of the matrix (gives a 3D tensor), and checks vs the analytical value. Note that 
 *  since the mesh is linear, it can only get the derivatives analytically correct for
 *  a linear function.
 */

void setMatField(apf::Mesh* mesh, apf::Field* nodal_fld)
{
  apf::MeshIterator* it=mesh->begin(0);
  apf::Field* coords_field = mesh->getCoordinateField();
  apf::MeshEntity* ent;
  apf::Vector3 coords;
  while((ent = mesh->iterate(it))) {
    apf::getVector(coords_field, ent, 0, coords);
    apf::Matrix3x3 F(1+5*coords[0], 0, 2*coords[1], 0, 1+coords[1], 0, 0, 0, 1+3*coords[2]);
    apf::setMatrix(nodal_fld, ent, 0, F);
  }
  mesh->end(it);
}

class MatrixDerivIntegrator : public apf::Integrator
{
  private:
    apf::MeshElement * mesh_elmt;
    apf::Vector3 global;
    // matrix field
    apf::Field * matrix_fld;
    // gradient of matrix field
    apf::Field * matrix_deriv_fld;
    apf::Vector<27> matrix_deriv;
    apf::Element * matrix_deriv_elmt;
  public:
    MatrixDerivIntegrator(int order, apf::Field* nodal_matrix_fld, apf::Field * matrix_deriv_fld) : Integrator(order), matrix_fld(nodal_matrix_fld), matrix_deriv_fld(matrix_deriv_fld) {}
    virtual void inElement(apf::MeshElement * elmt)
    {
      mesh_elmt = elmt;
      matrix_deriv_elmt = apf::createElement(matrix_fld, elmt);
    }
    virtual void atPoint(apf::Vector3 const& p, double, double)
    {
      apf::getMatrixGrad(matrix_deriv_elmt, p, matrix_deriv);
      apf::setComponents(matrix_deriv_fld, apf::getMeshEntity(mesh_elmt), ipnode, &matrix_deriv[0]);
    }
    virtual void outElement()
    {
      apf::destroyElement(matrix_deriv_elmt);
    }
};

class CheckMatrixDerivIntegrator : public apf::Integrator
{
  private:
    apf::MeshElement * mesh_elmt;
    apf::Vector3 global;
    apf::Field * matrix_deriv_fld;
  public:
    CheckMatrixDerivIntegrator(int order, apf::Field * matrix_deriv_fld) : Integrator(order), matrix_deriv_fld(matrix_deriv_fld) {}
    virtual void inElement(apf::MeshElement * elmt)
    {
      mesh_elmt = elmt;
    }
    virtual void atPoint(apf::Vector3 const& p, double, double)
    {
      apf::mapLocalToGlobal(mesh_elmt, p, global); 
      double dFdX[27] = {5, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 2, 0, 1, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 0, 0, 3};
      double values[27];
      apf::getComponents(matrix_deriv_fld, apf::getMeshEntity(mesh_elmt), ipnode, values);
      for(int i=0; i<27; ++i)
      {
        PCU_ALWAYS_ASSERT(fabs(values[i]-dFdX[i])<1E-10);
      }
    }
};

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    std::cerr<<"Usage: "<<argv[0]<<" model.dmg mesh.smb"<<std::endl;
    return 1;
  }
  MPI_Init(&argc, &argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  apf::Mesh2* mesh = apf::loadMdsMesh(argv[1], argv[2], &pcu_obj);
  int order=1;
  apf::Field* nodal_matrix_fld = apf::createLagrangeField(mesh, "matrix", apf::MATRIX, order);
  apf::Field* matrix_deriv = apf::createPackedField(mesh, "matrix_deriv", 27, apf::getIPShape(3,order));
  std::cout<<"loaded apf mesh created!"<<std::endl;
  apf::printStats(mesh);
  std::cout<<"Applying matrix field"<<std::endl;
  setMatField(mesh, nodal_matrix_fld);
  std::cout<<"Setting matrix derivatives"<<std::endl;
  apf::Integrator * set_matrix_deriv = new MatrixDerivIntegrator(order, nodal_matrix_fld, matrix_deriv);
  set_matrix_deriv->process(mesh);
  delete set_matrix_deriv;
  std::cout<<"Checking values"<<std::endl;
  apf::Integrator * check_matrix_deriv = new CheckMatrixDerivIntegrator(order, matrix_deriv);
  check_matrix_deriv->process(mesh);
  delete check_matrix_deriv;
  std::cout<<"Done"<<std::endl;
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  }
  MPI_Finalize();
  return 0;
}
