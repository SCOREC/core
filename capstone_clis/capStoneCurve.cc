#include <PCU.h>
#include <apfCAP.h>
#include <crv.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <lionPrint.h>
#include <memory>

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */


#include "CapstoneModule.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Analysis.h"
#include "CreateMG_Framework_Application.h"
#include "CreateMG_Framework_Attributes.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Geometry.h"
#include "CreateMG_Framework_Mesh.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

void safe_mkdir(const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("Err: could not create directory \"%s\"\n", path);
  }
}

void inspectFace(apf::Mesh2* m, apf::MeshEntity* e)
{
  apf::MeshEntity* down[12];
  int ndv = m->getDownward(e, 0, down);
  int nde = m->getDownward(e, 1, down);
  int ndf = m->getDownward(e, 2, down);
  printf("downward v/e/f %d/%d/%d\n", ndv, nde, ndf);
}

void inspectRegion(apf::Mesh2* m, apf::MeshEntity* e)
{
  apf::MeshEntity* down[12];
  int ndv = m->getDownward(e, 0, down);
  int nde = m->getDownward(e, 1, down);
  int ndf = m->getDownward(e, 2, down);
  int ndr = m->getDownward(e, 3, down);
  printf("downward v/e/f/r %d/%d/%d/%d\n", ndv, nde, ndf, ndr);
}

void writeVtk(CapstoneModule& cs, const std::string& vtkFileName)
{
  GeometryDatabaseInterface    *gdbi = cs.get_geometry();
  MeshDatabaseInterface        *mdbi = cs.get_mesh();
  AppContext		       *ctx = cs.get_context();

  // Get the VTK writer.
  Writer *vtkWriter = get_writer(ctx, "VTK File - Generic Writer");
  if (!vtkWriter)
	  error(HERE, ERR_GENERIC, "Could not find the VTK writer");

  // This IO Tool re-indexes the vertices for the elements classified in the gregions set and writes out the VTK file
  // the sets for edges and faces are empty so the mesh edge/faces classified over the egdes and faces wont get written explicitly

  /* std::vector<M_GTopo> vgregions; */
  /* vgregions.insert(vgregions.end(), gregions.begin(), gregions.end()); */
  /* vtkWriter->input().set_value("Regions", vgregions); */

  IdMapping idmapping;
  std::vector<M_MModel> mmodels;
  M_GModel gmodel;
  M_MModel mmodel;
  gdbi->get_current_model(gmodel);
  mdbi->get_current_model(mmodel);
  mmodels.clear();
  mmodels.push_back(mmodel);
  vtkWriter->write(ctx, gmodel, mmodels, vtkFileName.c_str(), idmapping);

  //cs.save_mesh_file(SystemTool::get_real_app_path().append("/booleaned.vtk"),ommodel);

}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));

  lion_set_verbosity(1);

  gmi_cap_start();
  gmi_register_cap();

  // load capstone mesh

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
  filenames.push_back(argv[1]);

  M_GModel gmodel = cs.load_files(filenames);

  M_MModel mmodel;
  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));

  for(std::size_t i = 0; i < mmodels.size(); ++i)
  {
      M_MModel ammodel = mmodels[i];
      std::size_t numregs = 0;
      MG_API_CALL(m, set_current_model(ammodel));
      MG_API_CALL(m, get_num_topos(TOPO_REGION, numregs));
      if(numregs > 0)
      {
	  mmodel = ammodel;
	  break;
      }
  }

  MG_API_CALL(m, set_adjacency_state(REGION2FACE|
				     REGION2EDGE|
				     REGION2VERTEX|
				     FACE2EDGE|
				     FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(m, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(m, compute_adjacency());

  /* writeVtk(cs, "before.vtk"); */
  apf::Mesh2* apfCapMesh = apf::createMesh(m, g, PCUObj.get());

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

  apf::writeVtkFiles("before", apfCapMesh);

  int order = 2;
  crv::BezierCurver bc(apfCapMesh, order, 0);
  bc.run();

  crv::writeCurvedVtuFiles(apfCapMesh, apf::Mesh::TRIANGLE, order + 2, "curved");
  crv::writeCurvedWireFrame(apfCapMesh, order + 10, "curved");


  M_MModel mmdl;
  m->get_current_model(mmdl);
  std::string info;
  m->print_info(mmdl, info);
  std::cout << info << std::endl;

  apf::writeVtkFiles("after", apfCapMesh);
  /* writeVtk(cs, "after.vtk"); */

  gmi_cap_stop();
  }
  MPI_Finalize();
}
