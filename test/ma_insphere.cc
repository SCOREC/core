#include <apfMDS.h>
#include <maMesh.h>
#include <gmi_null.h>
#include <PCU.h>
#include <apf.h>

int main(int argc, char** argv)
{
	MPI_Init(&argc,&argv);

	// Test determinant functions
	double input[4][4] = {
			{2, 5, 3, 5},
			{14, 9, 6, 7},
			{4, 9, 3, 2},
			{3, 7, 8, 6}
	};
  apf::Matrix<4,4> matrix;
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    matrix[i][j] = input[i][j];

	assert(apf::getDeterminant(apf::getMinor(matrix,0,0)) == 135);
	assert(apf::getDeterminant(apf::getMinor(matrix,1,0)) == 145);
	assert(apf::getDeterminant(apf::getMinor(matrix,2,0)) == 35);
	assert(apf::getDeterminant(apf::getMinor(matrix,3,0)) == -45);

	assert(apf::getDeterminant(matrix) == -1485);

	// Test insphere (create a mesh with one tet)
	PCU_Comm_Init();
	apf::Vector3 a(0, 0, 0);
	apf::Vector3 b(-6, 0, 0);
	apf::Vector3 c(0, -6, 0);
	apf::Vector3 d(0, 0, 12);

	gmi_register_null();
	gmi_model* model = gmi_load(".null");
	apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, 3, true);
	apf::ModelEntity* m = mesh->findModelEntity(0, 0);
	apf::MeshEntity* v[4];
	for (int i=0; i<4; i++) {
		v[i] = mesh->createVert(m);
	}
	mesh->setPoint(v[0], 0, a);
	mesh->setPoint(v[1], 0, b);
	mesh->setPoint(v[2], 0, c);
	mesh->setPoint(v[3], 0, d);

	apf::MeshEntity* e = apf::buildElement(mesh, m, apf::Mesh::TET, v);

	assert(ma::getInsphere(mesh, e) == 1.5);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

	PCU_Comm_Free();
	MPI_Finalize();
}
