#include <PCU.h>
#include <apf.h>
#include <apfShape.h>
#include <apfCAP.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <maSnap.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>

#include "CapstoneModule.h"
/* #include "CreateMG_Framework_Core.h" */
/* #include "CreateMG_Framework_Analysis.h" */
/* #include "CreateMG_Framework_Application.h" */
/* #include "CreateMG_Framework_Attributes.h" */
/* #include "CreateMG_Framework_Core.h" */
/* #include "CreateMG_Framework_Geometry.h" */
/* #include "CreateMG_Framework_Mesh.h" */


using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

class B737 : public ma::IsotropicFunction
{
  public:
    B737(ma::Mesh* m)
    {
      mesh = m;
    }
    virtual double getValue(ma::Entity* v)
    {
      // get current size
      apf::Up up;
      mesh->getUp(v, up);
      double sum = 0;
      double count = 0;
      for (int i = 0; i < up.n; i++) {
	ma::Entity* currentEdge = up.e[i];
	ma::Entity* vs[2];
	mesh->getDownward(currentEdge, 0, vs);
	ma::Vector diff = ma::getPosition(mesh, vs[0]) - ma::getPosition(mesh, vs[1]);
	sum += diff.getLength();
	count++;
      }
      double currentSize = sum/count;

      // compute the size multiplier
      double y = 1.;
      int tag  = mesh->getModelTag(mesh->toModel(v));
      int type = mesh->getModelType(mesh->toModel(v));
      if (type == 2 && (tag == 66 || tag == 86))
      	y = 3;
      if (type == 1 && tag == 178)
      	y = 3;
      // new size is multiplier * currentSize
      return y*currentSize;
    }
  private:
    ma::Mesh* mesh;
};


void writeCre(CapstoneModule& cs, const std::string& vtkFileName)
{
  GeometryDatabaseInterface    *gdbi = cs.get_geometry();
  MeshDatabaseInterface        *mdbi = cs.get_mesh();
  AppContext		       *ctx = cs.get_context();

  // Get the VTK writer.
  Writer *creWriter = get_writer(ctx, "Create Native Writer");
  if (!creWriter)
	  error(HERE, ERR_GENERIC, "Could not find the CRE writer");

  IdMapping idmapping;
  std::vector<M_MModel> mmodels;
  M_GModel gmodel;
  M_MModel mmodel;
  gdbi->get_current_model(gmodel);
  mdbi->get_current_model(mmodel);
  mmodels.clear();
  mmodels.push_back(mmodel);
  creWriter->write(ctx, gmodel, mmodels, vtkFileName.c_str(), idmapping);
}

void gradeSize(apf::Mesh2* m, apf::Field* f, double beta)
{
  std::queue<apf::MeshEntity*> q;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    apf::MeshEntity* vs[2];
    m->getDownward(e, 0, vs);
    double s0 = apf::getScalar(f, vs[0], 0);
    double s1 = apf::getScalar(f, vs[1], 0);
    double factor;
    if (s0 > s1) factor = s0/s1;
    else factor = s1/s0;
    if (factor > beta)
      q.push(e);
  }
  m->end(it);

  while(!q.empty()) {
    apf::MeshEntity* e = q.front();
    q.pop();

    apf::MeshEntity* vs[2];
    m->getDownward(e, 0, vs);
    double s0 = apf::getScalar(f, vs[0], 0);
    double s1 = apf::getScalar(f, vs[1], 0);

    apf::Up up;
    if (s1 > beta*s0) {
      s1 = beta*s0;
      apf::setScalar(f, vs[1], 0, s1);
      m->getUp(vs[1], up);
    }
    else if (s0 > beta*s1) {
      s0 = beta*s1;
      apf::setScalar(f, vs[0], 0, s0);
      m->getUp(vs[0], up);
    }
    else
      continue;

    for (int i = 0; i < up.n; i++) {
      q.push(up.e[i]);
    }
  }
}

static double triQuality(apf::Mesh2* m, apf::MeshEntity* tri)
{
  apf::Vector3 p0;
  apf::Vector3 p1;
  apf::Vector3 p2;

  apf::MeshEntity* vs[3];
  m->getDownward(tri, 0, vs);

  m->getPoint(vs[0], 0, p0);
  m->getPoint(vs[1], 0, p1);
  m->getPoint(vs[2], 0, p2);

  apf::Vector3 d1 = p1 - p0;
  apf::Vector3 d2 = p2 - p0;

  apf::Vector3 c = apf::cross(d1, d2);

  double areaSquare = c*c / 4.;

  apf::MeshEntity* es[3];
  m->getDownward(tri, 1, es);
  double sum = 0;
  for (int i = 0; i < 3; i++) {
    apf::MeshEntity* vs[2];
    m->getDownward(es[i], 0, vs);
    apf::Vector3 p0;
    apf::Vector3 p1;
    m->getPoint(vs[0], 0, p0);
    m->getPoint(vs[1], 0, p1);
    sum += (p1-p0)*(p1-p0);
  }
  return 48*areaSquare/sum/sum;
}


void stats(apf::Mesh2* m, apf::Field* f, double desiredQ, double &minQ, int &nQ, int &nMinL, int &nMaxL)
{
  apf::MeshEntity* e;
  apf::MeshIterator* it;

  nQ = nMinL = nMaxL = 0;

  minQ = 1000.0;
  it = m->begin(2);
  while ( (e = m->iterate(it)) ) {
    double q = triQuality(m, e);
    if (q < minQ)
      minQ = q;
    if (q < desiredQ)
      nQ++;
  }
  m->end(it);

  it = m->begin(1);
  while ( (e = m->iterate(it)) ) {
    apf::MeshEntity* vs[2];
    m->getDownward(e, 0, vs);
    apf::Vector3 p0;
    apf::Vector3 p1;
    m->getPoint(vs[0], 0, p0);
    m->getPoint(vs[1], 0, p1);
    double length = (p1-p0).getLength();
    double s0 = apf::getScalar(f, vs[0],0);
    double s1 = apf::getScalar(f, vs[1],0);
    double normalizedLength = 2*length/(s0+s1);
    if (normalizedLength < 0.5)
      nMinL++;
    if (normalizedLength > 1.5)
      nMaxL++;
  }
  m->end(it);
}

int main(int argc, char** argv)
{
  /* INIT CALLS */
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);

  if (argc < 2) {
    if (PCUObj.Self() == 0) {
      printf("USAGE: %s <in_cre_file.cre>\n", argv[0]);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  const char* inFileName    = argv[1];

  gmi_cap_start();
  gmi_register_cap();

  /* LOAD CAPSTONE MESH */
  // create an instance of the Capstone Module activating CREATE/CREATE/CREATE
  // for the Geometry/Mesh/Attribution databases
  const std::string gdbName("Geometry Database : SMLIB");// Switch Create with SMLIB for CAD
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("test", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();
  AppContext                    *c = cs.get_context();

  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);
  PCU_ALWAYS_ASSERT(c);

  v_string filenames;
  filenames.push_back(inFileName);

  M_GModel gmodel = cs.load_files(filenames);

  M_MModel mmodel;
  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
  PCU_ALWAYS_ASSERT(mmodels.size() == 1);
  MG_API_CALL(m, set_current_model(mmodels[0]));

  /* SET THE ADJACENCIES */
  MG_API_CALL(m, set_adjacency_state(REGION2FACE|
				     REGION2EDGE|
				     REGION2VERTEX|
				     FACE2EDGE|
				     FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(m, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(m, compute_adjacency());

  /* CONVERT THE MESH TO APF::MESH2 */
  apf::Mesh2* apfCapMesh = apf::createMesh(m, g, &PCUObj);

  /* ADD A TEST FIELD TO THE MESH TO DEMONSTRATE SOLUTION TRANSFER */
  apf::Field* tf  = apf::createFieldOn(apfCapMesh, "test_field", apf::VECTOR);
  apf::MeshEntity* ent;
  apf::MeshIterator* it = apfCapMesh->begin(0);
  while( (ent = apfCapMesh->iterate(it)) ) {
    apf::Vector3 p = ma::getPosition(apfCapMesh, ent);
    double x = p[0];
    double y = p[1];
    double z = p[2];
    apf::Vector3 s(y, z, 2*x);
    apf::setVector(tf, ent, 0, s);
  }
  apfCapMesh->end(it);

  /* WRITE THE BEFORE ADAPT MESH TO VTK USING APF VTK WRITER */
  apf::writeVtkFiles("before", apfCapMesh);


  /* SIZE SETUP AND CALL ADAPT */

  // setting the isotropic sizes
  ma::IsotropicFunction* sf = new B737(apfCapMesh);
  apf::Field* adaptSize  = apf::createFieldOn(apfCapMesh, "adapt_size", apf::SCALAR);
  apf::MeshEntity* v;
  it = apfCapMesh->begin(0);
  while( (v = apfCapMesh->iterate(it)) ) {
    apf::setScalar(adaptSize, v, 0, sf->getValue(v));
  }
  apfCapMesh->end(it);
  apf::writeVtkFiles("size_before_grade", apfCapMesh);

  // grading the size field
  gradeSize(apfCapMesh, adaptSize, 1.4);
  apf::writeVtkFiles("size_after_grade", apfCapMesh);

  // adapt setup
  ma::Input* in;
  in = ma::makeAdvanced(ma::configure(apfCapMesh, adaptSize));
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->shouldForceAdaptation = true;
  if (apfCapMesh->getDimension() == 2)
    in->goodQuality = 0.04;
  else // 3D meshes
    in->goodQuality = 0.027; // this is the mean-ratio cubed
  in->maximumIterations = 10;

  ma::adaptVerbose(in, false);


  /* CREATE A CONSTANT FIELD ON THE MESH WITH QUALITIES */
  apf::Field* qual  = apf::createField(apfCapMesh, "quality", apf::SCALAR, apf::getConstant(2));
  apf::MeshEntity* f;
  it = apfCapMesh->begin(2);
  while ( (f = apfCapMesh->iterate(it)) ) {
    double q = triQuality(apfCapMesh, f);
    apf::setScalar(qual, f, 0, q);
  }
  apfCapMesh->end(it);


  /* WRITE THE AFTER ADAPT MESH TO VTK USING APF VTK WRITER */
  apf::writeVtkFiles("after", apfCapMesh);


  /* PRINTING OUT SOME STATES */
  size_t vCount = apfCapMesh->count(0);
  size_t eCount = apfCapMesh->count(1);
  size_t fCount = apfCapMesh->count(2);

  double minQ;
  int nQ, nMinL, nMaxL;


  stats(apfCapMesh, adaptSize, 0.2*0.2, minQ, nQ, nMinL, nMaxL);

  printf("vCount is %d\n", (int) vCount);
  printf("eCount is %d\n", (int) eCount);
  printf("fCount is %d\n", (int) fCount);
  printf("minQ   is %f\n", minQ);
  printf("nQ     is %d\n", nQ);
  printf("nMinL  is %d\n", nMinL);
  printf("nMaxL  is %d\n", nMaxL);

  /* PRINT ADAPTED MESH INFO */
  M_MModel mmdl;
  m->get_current_model(mmdl);
  std::string info;
  m->print_info(mmdl, info);

  /* WRITE THE ADAPTED MESH IN NATIVE CREATE FORMAT */
  writeCre(cs, "after.cre");

  /* CLEAN UPS */
  delete sf;

  /* EXIT CALLS */
  gmi_cap_stop();
  }
  MPI_Finalize();
}
