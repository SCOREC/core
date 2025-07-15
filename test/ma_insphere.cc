#include <apfMDS.h>
#include <maMesh.h>
#include <gmi_null.h>
#include <lionPrint.h>
#include <apf.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
	pcu::Init(&argc,&argv);

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

	PCU_ALWAYS_ASSERT(apf::getDeterminant(apf::getMinor(matrix,0,0)) == 135);
	PCU_ALWAYS_ASSERT(apf::getDeterminant(apf::getMinor(matrix,1,0)) == 145);
	PCU_ALWAYS_ASSERT(apf::getDeterminant(apf::getMinor(matrix,2,0)) == 35);
	PCU_ALWAYS_ASSERT(apf::getDeterminant(apf::getMinor(matrix,3,0)) == -45);

	PCU_ALWAYS_ASSERT(apf::getDeterminant(matrix) == -1485);

	// Test insphere (create a mesh with one tet)
	{
	pcu::PCU PCUObj;
  lion_set_verbosity(1);
	apf::Vector3 a(0, 0, 0);
	apf::Vector3 b(-6, 0, 0);
	apf::Vector3 c(0, -6, 0);
	apf::Vector3 d(0, 0, 12);

	gmi_register_null();
	gmi_model* model = gmi_load(".null");
	apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, 3, true, &PCUObj);
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

	PCU_ALWAYS_ASSERT(ma::getInsphere(mesh, e) == 1.5);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

	}
	pcu::Finalize();
}
