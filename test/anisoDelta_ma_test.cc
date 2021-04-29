#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include "pumi.h"
#include <apfMatrix.h>

#include <stdlib.h>
#include "maStats.h"
#include <iostream>
#include <fstream>
#include <sstream>

//TODO commented blocks are for setting the input size and frames
/*
void getValue(ma::Mesh *m, ma::Entity* v, ma::Matrix& R, ma::Vector& H) {
  auto targetMetric = m->findField("target_metric");
  PCU_ALWAYS_ASSERT(targetMetric);
  int nComps = targetMetric->countComponents();
  auto vals = new double[nComps];
  apf::getComponents(targetMetric, v, 0, vals);
  auto scale = 1e4;
  vals[0] = vals[0]/scale;
  vals[1] = vals[1]/scale;
  vals[2] = vals[2]/scale;
  vals[3] = vals[3]/scale;
  vals[4] = vals[4]/scale;
  vals[5] = vals[5]/scale;
  ma::Matrix M(vals[0], vals[3], vals[5],
               vals[3], vals[1], vals[4],
               vals[5], vals[4], vals[2]);

  double eigenValues[3];
  ma::Vector eigenVectors[3];
  apf::eigen(M, eigenVectors, eigenValues);
  ma::Matrix RT(eigenVectors); // eigen vectors are stored in the rows of RT
  RT[0] = RT[0].normalize();
  RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
  RT[1] = RT[1].normalize();
  RT[2] = apf::cross(RT[0],RT[1]);
  R = apf::transpose(RT);

  double h[3];
  for (int i = 0; i < 3; ++i)
    h[i] = std::sqrt(1.0/(eigenValues[i]*scale));
  H = ma::Vector(h);
  delete [] vals;
}
*/

class AnIso : public ma::AnisotropicFunction {
public:
  AnIso(ma::Mesh* m) {
    mesh = m;
    ma::getBoundingBox(m,lower,upper);
    targetMetric = m->findField("target_metric");
    PCU_ALWAYS_ASSERT(targetMetric);
    int nComps = targetMetric->countComponents();
    vals = new double[nComps];
  }
  ~AnIso() {
    delete [] vals;
  }
  virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H) {
    apf::getComponents(targetMetric, v, 0, vals);
    auto scale = 1e4;
    vals[0] = vals[0]/scale;
    vals[1] = vals[1]/scale;
    vals[2] = vals[2]/scale;
    vals[3] = vals[3]/scale;
    vals[4] = vals[4]/scale;
    vals[5] = vals[5]/scale;
    ma::Matrix M(vals[0], vals[3], vals[5],
                 vals[3], vals[1], vals[4],
	         vals[5], vals[4], vals[2]);

    double eigenValues[3];
    ma::Vector eigenVectors[3];
    apf::eigen(M, eigenVectors, eigenValues);
    ma::Matrix RT(eigenVectors); // eigen vectors are stored in the rows of RT
    RT[0] = RT[0].normalize();
    RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
    RT[1] = RT[1].normalize();
    RT[2] = apf::cross(RT[0],RT[1]);
    R = apf::transpose(RT);

    double h[3];
    for (int i = 0; i < 3; ++i)
      h[i] = std::sqrt(1.0/(eigenValues[i]*scale));
    H = ma::Vector(h);
  }
private:
  ma::Mesh* mesh;
  ma::Vector lower;
  ma::Vector upper;
  apf::Field* targetMetric;
  double* vals;
};

int main(int argc, char** argv) {
  PCU_ALWAYS_ASSERT(argc == 4);
  const char *modelFile = argv[1];
  const char *meshFile = argv[2];
  bool logInterpolation = true;
  MPI_Init(&argc, &argv);
  PCU_Comm_Init ();
  lion_set_verbosity (1);
  gmi_register_mesh ();
  gmi_register_null ();
  ma::Mesh* m = apf::loadMdsMesh (modelFile, meshFile);
  auto targetMetric = m->findField ("target_metric");
  PCU_ALWAYS_ASSERT (targetMetric);
  m->verify();
  apf::writeVtkFiles ("anisoDelta_before",m);

/*
  apf::MeshEntity* ent;
  apf::MeshIterator* it; 

  apf::Field* sizes = apf::createField(m, "sizes", apf::VECTOR, apf::getLagrange(1));
  apf::Field* frames = apf::createField(m, "frames", apf::MATRIX, apf::getLagrange(1));
  it = m->begin(0);
  while( (ent = m->iterate(it)) ) {
    apf::Vector3 sz; // use what your have to get this H
    apf::Matrix3x3 frm; // use what you have to get this R
    getValue(ent, frm, sz);
    apf::setVector(sizes, ent, 0, sz);
    apf::setMatrix(frames, ent, 0, frm);
  }
  m->end(it);

  ma::Input *in = ma::configure(m, sizes, frames, 0, true);
*/

  AnIso sf (m);
  // the metrics should be given as input using sizes and frames same as
  // test/highOrderSizeFields.cc
  ma::Input* in = ma::configure (m, &sf, 0, logInterpolation);
  in-> shouldRunPreZoltan = true;
  in-> shouldRunMidParma = true;
  in-> shouldRunPostParma = true;
  in-> shouldRefineLayer = true;
  in-> goodQuality = 0.027;

  in->maximumIterations = 10;

  ma::adaptVerbose (in);
  m-> verify ();

  std::vector<double> lengths;
  std::vector<double> qualities;
  in = ma::configure (m, &sf, 0, logInterpolation);
  ma::stats(m, in->sizeField, lengths, qualities, true);// we need to define
  // the size fields in order to get the lengths and qualities

  std::ostringstream len_fileNameStream("lengths_");
  len_fileNameStream << PCU_Comm_Self() << ".txt";
  std::string len_fileName = len_fileNameStream.str();
  std::ofstream len_file;
  len_file.open (len_fileName.c_str());
  for(std::size_t i = 0; i < lengths.size(); ++i) {
    len_file << lengths[i] << "\n";
  }
  len_file.close();

  std::ostringstream qua_fileNameStream("qualities_");
  qua_fileNameStream << PCU_Comm_Self() << ".txt";
  std::string qua_fileName = qua_fileNameStream.str();
  std::ofstream qua_file;
  qua_file.open (qua_fileName.c_str());
  for(std::size_t i = 0; i < qualities.size(); ++i) {
    qua_file << qualities[i] << "\n";
  }
  qua_file.close();

  apf::writeVtkFiles ("anisoDelta_after",m);

  m->destroyNative ();
  apf::destroyMesh (m);
  PCU_Comm_Free ();
  MPI_Finalize ();
}
