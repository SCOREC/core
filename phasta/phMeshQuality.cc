#include <phastaChef.h>
#include <gmi_mesh.h>

#include "apfSIM.h"
#include "gmi_sim.h"
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfMatrix.h>
#include <maSize.h>
#include <maShape.h>
#include <SimUtil.h>
#include <SimModel.h>
#include <SimPartitionedMesh.h>
#include <SimMeshTools.h>

#include <pcu_util.h>
#include <cstdio>

namespace ph {

void attachSIMSizeField(apf::Mesh2* m, apf::Field* sf_mag, apf::Field* sf_dir) {
// loop over all vertices
  double size[1];
  double anisosize[3][3];
  apf::MeshEntity* v;
  apf::MeshIterator* vit = m->begin(0);
  while ((v = m->iterate(vit))) {
// get sim size field
// this is for simmetrix mesh, should be generalized
    pVertex meshVertex = reinterpret_cast<pVertex>(v);
    int sztype = V_size(meshVertex, size, anisosize);
    PCU_ALWAYS_ASSERT(sztype == 1 || sztype == 2);

// consider iso size field as anisotropic size field
    if (sztype == 1) {
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          if (i==j)
            anisosize[i][j] = size[0];
          else
            anisosize[i][j] = 0.0;
        }
      }
    }

// transfer to apf field sf_mag and sf_dir
/* note that the frame in Simmetrix is stored by row
   while in PUMI, it is stored by column */
    apf::Vector3 v_mag;
    v_mag[0] = sqrt(anisosize[0][0]*anisosize[0][0]
                  + anisosize[0][1]*anisosize[0][1]
                  + anisosize[0][2]*anisosize[0][2]);
    v_mag[1] = sqrt(anisosize[1][0]*anisosize[1][0]
                  + anisosize[1][1]*anisosize[1][1]
                  + anisosize[1][2]*anisosize[1][2]);
    v_mag[2] = sqrt(anisosize[2][0]*anisosize[2][0]
                  + anisosize[2][1]*anisosize[2][1]
                  + anisosize[2][2]*anisosize[2][2]);
    apf::Matrix3x3 v_dir;
    v_dir = apf::Matrix3x3(
             anisosize[0][0]/v_mag[0], anisosize[1][0]/v_mag[1], anisosize[2][0]/v_mag[2],
             anisosize[0][1]/v_mag[0], anisosize[1][1]/v_mag[1], anisosize[2][1]/v_mag[2],
             anisosize[0][2]/v_mag[0], anisosize[1][2]/v_mag[1], anisosize[2][2]/v_mag[2]);
    apf::setVector(sf_mag, v, 0, v_mag);
    apf::setMatrix(sf_dir, v, 0, v_dir);
  }
  m->end(vit);
}

} // end namespace ph

#ifdef __cplusplus
extern"C"{
#endif

/* input x1,x2,x3 are current coordinates of mesh in phasta
   output minq is the minimum quality of mesh */
void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double& minq) {
// loop over all vertices
  int counter = 0;
  apf::MeshEntity* v;
  apf::MeshIterator* vit = m->begin(0);
  while ((v = m->iterate(vit))) {
// update the coordinates of current mesh
    apf::Vector3 p;
    p[0] = x1[counter];
    p[1] = x2[counter];
    p[2] = x3[counter];
    m->setPoint(v, 0, p);
    counter++;
  }
  m->end(vit);
  PCU_ALWAYS_ASSERT(counter==numnp);

// transfer to apf sizefield
  if(m->findField("sizes")) apf::destroyField(m->findField("sizes"));
  if(m->findField("frames")) apf::destroyField(m->findField("frames"));
// this is for simmetrix mesh, should be generalized
  apf::Field* sf_mag = apf::createSIMFieldOn(m, "sizes", apf::VECTOR);
  apf::Field* sf_dir = apf::createSIMFieldOn(m, "frames", apf::MATRIX);
  ph::attachSIMSizeField(m, sf_mag, sf_dir);
  ma::SizeField* sf = ma::makeSizeField(m, sf_mag, sf_dir);

// measure
  minq = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* eit = m->begin(m->getDimension());
  while( (e = m->iterate(eit)) ) {
    if (! m->isOwned(e))
      continue;
    double q = ma::measureElementQuality(m, sf, e);
    if (q < minq)
      minq = q;
  }
  m->end(eit);
}

#ifdef __cplusplus
}
#endif

