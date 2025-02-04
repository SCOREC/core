#include <apf.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>
#include <gmi_mesh.h>
#include <ma.h>
#include <maShape.h>
#include <crv.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <SimUtil.h>
#include <MeshSim.h>
#include <gmi_sim.h>
#include <SimModel.h>
#endif
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>


void computeSizesFrames(
    apf::Mesh2* m,
    const apf::Vector3 &p,
    apf::Vector3 &sz,
    apf::Matrix3x3 &frm);

void testAdapt(
    const char* model,
    const char* mesh,
    int order,
    pcu::PCU *PCUObj);

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);

#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();

  if (argc != 3) {
    if(0==PCUObj.Self())
      std::cerr <<"usage: " << argv[0]
      	<< " <model.smd> <mesh.smb> <sizefield>\n";
    return EXIT_FAILURE;
  }

  for (int order = 3; order < 5; order++) {
    if(0==PCUObj.Self())
      lion_oprint(1, "Testing aniso adapt w/ sizefield order %d\n", order);
    testAdapt(argv[1], argv[2], order, &PCUObj);
  }

  }
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  MPI_Finalize();
}

void computeSizesFrames(
    apf::Mesh2* m,
    const apf::Vector3 &p,
    apf::Vector3 &sz,
    apf::Matrix3x3 &frm)
{
  apf::Vector3 llb;
  apf::Vector3 umb;
  ma::getBoundingBox(m, llb, umb);
  const double pi = 3.14152;
  apf::Vector3 diff = umb - llb;
  diff[1] = 0.0;
  // total length in the direction of the wave
  double length = diff.getLength();
  // normal in the direction of the wave
  apf::Vector3 normal = diff / length;
  apf::Vector3 tangent(-normal[2], 0., normal[0]);

  double maxSz = length / 12.;
  double minSz = length / 200.;
  double wavelength = length / 4.;


  double alpha = (maxSz - minSz)/2.;
  double beta =  (maxSz + minSz)/2.;

  double s = normal * p;

  sz = apf::Vector3(
      length/4.,
      length/4.,
      beta + alpha * sin(2*pi*s/wavelength)
      );

  frm[0] = normal;
  frm[1] = apf::Vector3(0., 1., 0.);
  frm[2] = tangent;
}

void testAdapt(
    const char* model,
    const char* mesh,
    int order,
    pcu::PCU *PCUObj)
{
  gmi_model* g = gmi_load(model);
  apf::Mesh2* m = apf::loadMdsMesh(g,mesh,PCUObj);
  m->verify();

  const ma::Input* in = ma::configureUniformRefine(m, 1, 0);
  ma::adapt(in);

  apf::FieldShape* fs = apf::getH1Shape(order);
  apf::Field* sizes = apf::createField(m, "sizes", apf::VECTOR, fs);
  apf::Field* frames = apf::createField(m, "frames", apf::MATRIX, fs);

  int dim = m->getDimension();
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!fs->countNodesOn(apf::Mesh::simplexTypes[d]))
      continue;
    it = m->begin(d);
    while( (ent = m->iterate(it)) ) {
      int type = m->getType(ent);
      int non = fs->countNodesOn(type);
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      for (int i = 0; i < non; i++) {
	apf::Vector3 xi, p, value;
	fs->getNodeXi(type, i, xi);
	apf::mapLocalToGlobal(me, xi, p);
	apf::Vector3 sz;
	apf::Matrix3x3 frm;
	computeSizesFrames(m, p, sz, frm);
	apf::setVector(sizes, ent, i, sz);
	apf::setMatrix(frames, ent, i, frm);
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }

  ma::Input* inAdv = ma::makeAdvanced(ma::configure(m, sizes, frames, 0, true));
  inAdv->shouldFixShape = true;
  inAdv->maximumIterations = 10;
  inAdv->shouldForceAdaptation = true;

  std::stringstream ss;
  ss << "before_adapt_with_ho_sizefield_order_" << order;
  apf::writeVtkFiles(ss.str().c_str(), m);
  ss.str("");
  ma::adaptVerbose(inAdv);
  ss << "after_adapt_with_ho_sizefield_order_" << order;
  apf::writeVtkFiles(ss.str().c_str(), m);

  m->removeField(sizes);
  m->removeField(frames);

  apf::destroyField(sizes);
  apf::destroyField(frames);

  m->destroyNative();
  apf::destroyMesh(m);
}
