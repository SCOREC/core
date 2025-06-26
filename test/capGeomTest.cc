#include <iostream>
#include <iomanip>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <lionPrint.h>
#include <pcu_util.h>

char const* const typeName[4] =
{"vertex",
 "  edge",
 "  face",
 "region"};

void printInfo(gmi_model* model, int dim);

void visualizeFace(gmi_model* model, gmi_ent* entity, int n, int m, const char* fileName, pcu::PCU *PCUObj);
void visualizeEdge(gmi_model* model, gmi_ent* entity, int n, const char* fileName, pcu::PCU *PCUObj);
void visualizeEdges(gmi_model* model, int n, const char* fileName, pcu::PCU *PCUObj);

int main(int argc, char** argv)
{
  pcu::Init(&argc, &argv);
  { // pcu object scope
  pcu::PCU PCUObj;

  PCU_ALWAYS_ASSERT(argc == 2);

  gmi_register_mesh();
  gmi_register_null();
  gmi_cap_start();
  gmi_register_cap();

  gmi_model* model = gmi_load(argv[1]);
  apf::Mesh2* mesh0 = apf::createCapMesh(model, &PCUObj);
  apf::disownCapModel(mesh0);
  apf::writeVtkFiles("mesh_no_param", mesh0);

  printInfo(model, 0);
  printf("\n");
  printInfo(model, 1);
  printf("\n");
  printInfo(model, 2);
  printf("\n");
  printInfo(model, 3);

    // faces in separate meshes named by tags
  gmi_ent* ge;
  gmi_iter* gi;

  gi = gmi_begin(model, 2);
  while( (ge = gmi_next(model, gi)) ){
    std::stringstream name_str;
    name_str << "face_" << gmi_tag(model, ge) << "_mesh";
    visualizeFace(model, ge, 100, 100, name_str.str().c_str(), &PCUObj);
  }
  gmi_end(model, gi);

  printf("------------------------\n");
  printf("creating mesh with param field\n");


  apf::Mesh2* mesh = apf::createCapMesh(model, &PCUObj);
  apf::Field* pf  = apf::createFieldOn(mesh, "param_field", apf::VECTOR);
  apf::Field* idf  = apf::createFieldOn(mesh, "id", apf::SCALAR);
  apf::MeshEntity* e;
  apf::MeshIterator* it = mesh->begin(0);
  while ( (e = mesh->iterate(it)) ) {
    apf::Vector3 param;
    mesh->getParam(e, param);
    double tag = (double) mesh->getModelTag(mesh->toModel(e));
    apf::setVector(pf, e, 0, param);
    apf::setScalar(idf, e, 0, tag);
  }
  mesh->end(it);

  apf::writeVtkFiles("mesh_with_param", mesh);

  gmi_cap_stop();

  } // pcu object scope
  pcu::Finalize();
}

void printInfo(gmi_model* model, int dim)
{
  int type = dim;
  std::stringstream ss;
  int isP0 = 0;
  int isP1 = 0;
  double range0[2];
  double range1[2];
  printf("print info for model ents with dimension %d\n", dim);
  printf("-------------------------------------------\n");

  gmi_ent* ge;
  gmi_iter* gi = gmi_begin(model, type);
  while( (ge = gmi_next(model, gi)) ){
    ss.str("");
    ss << "ent: "
       << "tag " << std::setw(4) << gmi_tag(model, ge) << " "
       << "type " << typeName[type] << " ";
    switch (type){
      case 0:
      	isP0 = 0;
      	isP1 = 0;
      	range0[0] = -1.;
      	range0[1] = -1.;
      	range1[0] = -1.;
      	range1[1] = -1.;
	break;
      case 1:
      	isP0 = gmi_periodic(model, ge, 0) ? 1:0;
      	isP1 = 0;
	gmi_range(model, ge, 0, range0);
      	range1[0] = -1.;
      	range1[1] = -1.;
	break;
      case 2:
      	isP0 = gmi_periodic(model, ge, 0) ? 1:0 ;
      	isP1 = gmi_periodic(model, ge, 1) ? 1:0 ;
	gmi_range(model, ge, 0, range0);
	gmi_range(model, ge, 1, range1);
	break;
      case 3:
      	isP0 = 0;
      	isP1 = 0;
      	range0[0] = -1.;
      	range0[1] = -1.;
      	range1[0] = -1.;
      	range1[1] = -1.;
	break;
    }
    ss << "isP " << "(" << isP0 << ", " << isP1 << ") "
       << "range " << "(" << range0[0] << ", " << range0[1] << ")--(" << range1[0] << ", " << range1[1] << ")";
    if (type == 1) {
      double param[2] = {0., 0.};
      param[0] = range0[0];
      double posi[3];
      gmi_eval(model, ge, param, &posi[0]);
      apf::Vector3 posiv(posi[0], posi[1], posi[2]);
      ss << " ++ " << posiv;
      param[0] = range0[1];
      gmi_eval(model, ge, param, &posi[0]);
      posiv = apf::Vector3(posi[0], posi[1], posi[2]);
      ss << " -- " << posiv;
    }
    ss << std::endl;
    printf("%s", ss.str().c_str());
  }
  gmi_end(model, gi); // end the iterator
  printf("-------------------------------------------\n");
}

void visualizeFace(gmi_model* model, gmi_ent* entity, int n, int m, const char* fileName, pcu::PCU *PCUObj)
{
  // assert type is 2
  // get the range
  double u_range[2];
  double v_range[2];
  gmi_range(model, entity, 0, &u_range[0]);
  gmi_range(model, entity, 1, &v_range[0]);
  // update the v_range by tol
  /* double tol = (v_range[1] - v_range[0]) / (m-1); */
  /* v_range[0] += tol; */
  /* v_range[1] -= tol; */
  double du = (u_range[1] - u_range[0]) / (n-1);
  double dv = (v_range[1] - v_range[0]) / (m-1);

  // make the array of vertex coordinates in the physical space
  std::vector<apf::Vector3> ps;
  std::vector<apf::Vector3> uvs;
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      double params[2];
      params[0] = u_range[0] + i * du;
      params[1] = v_range[0] + j * dv;
      double position[3];
      gmi_eval(model, entity, &params[0], &position[0]);
      apf::Vector3 p(position[0], position[1], position[2]);
      ps.push_back(p);
      apf::Vector3 uv(params[0], params[1], 0.);
      uvs.push_back(uv);
    }
  }

  // make the vertexes and set the coordinates using the array
  std::vector<apf::MeshEntity*> vs;
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false, PCUObj);
  for (size_t i = 0; i < ps.size(); i++) {
    apf::MeshEntity* vert = mesh->createVert(0);
    mesh->setPoint(vert, 0, ps[i]);
    vs.push_back(vert);
  }

  assert(vs.size() == ps.size());


  apf::MeshEntity* v[3];
  // make the lower/upper t elems
  for (int i = 0; i < n-1; i++) {
    for (int j = 0; j < m-1; j++) {
      // upper triangle
      v[0] = vs[(i + 0) + n * (j + 0)];
      v[1] = vs[(i + 0) + n * (j + 1)];
      v[2] = vs[(i + 1) + n * (j + 0)];
      apf::buildElement(mesh, 0, apf::Mesh::TRIANGLE, v);
      // upper triangle
      v[0] = vs[(i + 0) + n * (j + 1)];
      v[1] = vs[(i + 1) + n * (j + 1)];
      v[2] = vs[(i + 1) + n * (j + 0)];
      apf::buildElement(mesh, 0, apf::Mesh::TRIANGLE, v);
    }
  }

  apf::deriveMdsModel(mesh);
  mesh->acceptChanges();
  mesh->verify();
  apf::printStats(mesh);

  // parametric coordinates
  apf::Field* f = apf::createFieldOn(mesh, "param_coords", apf::VECTOR);
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  it = mesh->begin(0);
  int count = 0;
  while ( (e = mesh->iterate(it)) ) {
    apf::setVector(f, e, 0, uvs[count]);
    count++;
  }

  // tangent vectors
  apf::Field* ut = apf::createFieldOn(mesh, "u_tangent", apf::VECTOR);
  apf::Field* vt = apf::createFieldOn(mesh, "v_tangent", apf::VECTOR);
  apf::Field* nv = apf::createFieldOn(mesh, "normal", apf::VECTOR);

  it = mesh->begin(0);
  count = 0;
  while ( (e = mesh->iterate(it)) ) {
    double uTangent[3];
    double vTangent[3];
    double normal[3];
    double param[2] = {uvs[count].x(), uvs[count].y()};
    gmi_first_derivative(model, entity, param, uTangent, vTangent);
    gmi_normal(model, entity, param, normal);
    apf::setVector(ut, e, 0, apf::Vector3(uTangent[0], uTangent[1], uTangent[2]));
    apf::setVector(vt, e, 0, apf::Vector3(vTangent[0], vTangent[1], vTangent[2]));
    apf::setVector(nv, e, 0, apf::Vector3(normal[0], normal[1], normal[2]));
    count++;
  }

  apf::writeVtkFiles(fileName,mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);


}

void visualizeEdge(gmi_model* model, gmi_ent* entity, int n, const char* fileName, pcu::PCU *PCUObj)
{
  // assert type is 1
  // get the range
  double u_range[2];
  gmi_range(model, entity, 0, &u_range[0]);
  double du = (u_range[1] - u_range[0]) / (n-1);

  // make the array of vertex coordinates in the physical space
  std::vector<apf::Vector3> ps;
  std::vector<apf::Vector3> us;
  for (int i = 0; i < n; i++) {
    double params[2];
    params[0] = u_range[0] + i * du;
    double position[3];
    gmi_eval(model, entity, &params[0], &position[0]);
    apf::Vector3 p(position[0], position[1], position[2]);
    ps.push_back(p);
    apf::Vector3 uv(params[0], params[1], 0.);
    us.push_back(uv);
  }

  // make the vertexes and set the coordinates using the array
  std::vector<apf::MeshEntity*> vs;
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false, PCUObj);
  for (size_t i = 0; i < ps.size(); i++) {
    apf::MeshEntity* vert = mesh->createVert(0);
    mesh->setPoint(vert, 0, ps[i]);
    vs.push_back(vert);
  }

  assert(vs.size() == ps.size());


  apf::MeshEntity* v[2];
  // make the lower/upper t elems
  for (int i = 0; i < n-1; i++) {
    v[0] = vs[i];
    v[1] = vs[i+1];
    apf::buildElement(mesh, 0, apf::Mesh::EDGE, v);
  }

  apf::deriveMdsModel(mesh);
  mesh->acceptChanges();
  mesh->verify();
  apf::printStats(mesh);

  apf::Field* f = apf::createFieldOn(mesh, "param_coords", apf::VECTOR);
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  it = mesh->begin(0);
  int count = 0;
  while ( (e = mesh->iterate(it)) ) {
    apf::setVector(f, e, 0, us[count]);
    count++;
  }

  apf::writeVtkFiles(fileName,mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
}

void visualizeEdges(gmi_model* model, int n, const char* fileName, pcu::PCU *PCUObj)
{
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false, PCUObj);
  gmi_ent* entity;
  gmi_iter* gi = gmi_begin(model, 1);
  while( (entity = gmi_next(model, gi)) ){
    // assert type is 1
    // get the range
    double u_range[2];
    gmi_range(model, entity, 0, &u_range[0]);
    double du = (u_range[1] - u_range[0]) / (n-1);

    // make the array of vertex coordinates in the physical space
    std::vector<apf::Vector3> ps;
    std::vector<apf::Vector3> us;
    for (int i = 0; i < n; i++) {
      double params[2];
      params[0] = u_range[0] + i * du;
      double position[3];
      gmi_eval(model, entity, &params[0], &position[0]);
      apf::Vector3 p(position[0], position[1], position[2]);
      ps.push_back(p);
      apf::Vector3 uv(params[0], params[1], 0.);
      us.push_back(uv);
    }

    // make the vertexes and set the coordinates using the array
    std::vector<apf::MeshEntity*> vs;
    for (size_t i = 0; i < ps.size(); i++) {
      apf::MeshEntity* vert = mesh->createVert(0);
      mesh->setPoint(vert, 0, ps[i]);
      vs.push_back(vert);
    }

    assert(vs.size() == ps.size());

    apf::MeshEntity* v[2];
    // make the lower/upper t elems
    for (int i = 0; i < n-1; i++) {
      v[0] = vs[i];
      v[1] = vs[i+1];
      apf::buildElement(mesh, 0, apf::Mesh::EDGE, v);
    }
  }
  gmi_end(model, gi);


  apf::deriveMdsModel(mesh);
  mesh->acceptChanges();
  mesh->verify();
  apf::printStats(mesh);


  apf::writeVtkFiles(fileName,mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
}
