#include <phastaChef.h>
#include <gmi_mesh.h>

#include "apfSIM.h"
#include "gmi_sim.h"
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfMatrix.h>
#include <maSize.h>
#include <maShape.h>
//#include "apf_simConfig.h"

#ifdef HAVE_SIMMETRIX
#include <SimUtil.h>
#include <SimModel.h>
#include <SimPartitionedMesh.h>
#include <SimMeshTools.h>
#ifdef HAVE_SIMADVMESHING
  #include <SimAdvMeshing.h>
#endif
#endif

#include <PCU.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdio>

namespace ph {

static void initArray(double array[3][3]) {
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      array[i][j] = 0.0;
    }
  }
}

void attachSIMSizeField(apf::Mesh2* m, apf::Field* sf_mag, apf::Field* sf_dir) {
// loop over all vertices
  double size[1];
  double anisosize[3][3];
  apf::MeshEntity* v;
  apf::MeshIterator* vit = m->begin(0);
  while ((v = m->iterate(vit))) {
#ifdef HAVE_SIMMETRIX
// get sim size field
// this is for simmetrix mesh, should be generalized
    {
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      int sztype = V_size(meshVertex, size, anisosize);
      PCU_ALWAYS_ASSERT(sztype == 1 || sztype == 2);

      // consider iso size field as anisotropic size field
      if (sztype == 1) {
        initArray(anisosize);
        anisosize[0][0] = size[0];
        anisosize[1][1] = size[0];
        anisosize[2][2] = size[0];
      }
    }
#else
    {
      PCU_ALWAYS_ASSERT_VERBOSE(0,
         "please turn on Simmetrix and re-compile the code.\n");
    }
#endif

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

// the same as above but for isotropic size fields
void attachSIMSizeField(apf::Mesh2* m, apf::Field* sf_mag) {
// loop over all vertices
  double size[1];
  double anisosize[3][3];
  apf::MeshEntity* v;
  apf::MeshIterator* vit = m->begin(0);
  while ((v = m->iterate(vit))) {
#ifdef HAVE_SIMMETRIX
// get sim size field
// this is for simmetrix mesh, should be generalized
    {
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      int sztype = V_size(meshVertex, size, anisosize);
      PCU_ALWAYS_ASSERT(sztype == 1 || sztype == 2);

      // consider iso size field as anisotropic size field
      if (sztype == 1) {
        initArray(anisosize);
        anisosize[0][0] = size[0];
        anisosize[1][1] = size[0];
        anisosize[2][2] = size[0];
      }
    }
#else
    {
      PCU_ALWAYS_ASSERT_VERBOSE(0,
         "please turn on Simmetrix and re-compile the code.\n");
    }
#endif

// transfer to apf field sf_mag
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
    apf::setVector(sf_mag, v, 0, v_mag);
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
                        double& minvq, double& minfq) {
  double** op = (double**)malloc(sizeof(double*) * numnp);
  for (int i = 0; i < numnp; i++)
    op[i] = (double*)malloc(sizeof(double) * 3);
// loop over all vertices
  int counter = 0;
  apf::Vector3 p;
  apf::MeshEntity* v;
  apf::MeshIterator* vit = m->begin(0);
  while ((v = m->iterate(vit))) {
// get original mesh coordinates
    m->getPoint(v, 0, p);
    op[counter][0] = p[0];
    op[counter][1] = p[1];
    op[counter][2] = p[2];
// update the coordinates of current mesh
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
  apf::Field* sf_mag = apf::createSIMFieldOn(m, "sizes", apf::VECTOR);
  apf::Field* sf_dir = apf::createSIMFieldOn(m, "frames", apf::MATRIX);
// this is for simmetrix mesh, should be generalized
  PCU_ALWAYS_ASSERT_VERBOSE(in.simmetrixMesh, "core_measure_mesh only supports Simmetrix mesh!\n");
  ph::attachSIMSizeField(m, sf_mag, sf_dir);
  ma::SizeField* sf = ma::makeSizeField(m, sf_mag, sf_dir);

// measure mesh region
  minvq = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* eit = m->begin(m->getDimension());
  while( (e = m->iterate(eit)) ) {
    if (! m->isOwned(e))
      continue;
    /* not support other mesh region type */
    if (m->getType(e) != apf::Mesh::TET)
      continue;
#ifdef HAVE_SIMADVMESHING
    if (EN_isBLEntity(reinterpret_cast<pEntity>(e)))
      continue;
#endif
    double vq = ma::measureElementQuality(m, sf, e); // cubic mean ratio
    vq = cbrt(vq); // mean ratio
    if (vq < minvq)
      minvq = vq;
  }
  m->end(eit);

// measure mesh face in layered mesh
  minfq = 1.0;
#ifdef HAVE_SIMADVMESHING
  if (in.simmetrixMesh == 1) {
// get simmetrix mesh
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh pmesh = apf_msim->getMesh();
    pMesh mesh = PM_mesh(pmesh,0);

// get simmetrix model
    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);

// measure triangles
    pGFace gFace;
    pFace meshFace;
    pEntity seedRegion;
    pPList growthRegion = PList_new();
    pPList growthLayerFace = PList_new();
    GFIter gFIter = GM_faceIter(model);
    while((gFace = GFIter_next(gFIter))){
      FIter fIter = M_classifiedFaceIter(mesh, gFace, 1);
      while((meshFace = FIter_next(fIter))){
        if(BL_isBaseEntity(meshFace,gFace) == 1){
          for (int fromSide = 0; fromSide < 2; fromSide++) {
            int hasSeed = BL_stackSeedEntity(meshFace,gFace,fromSide,NULL,&seedRegion);
            PCU_ALWAYS_ASSERT_VERBOSE(hasSeed >= 0, "BL blending is not supported!\n");
            if (hasSeed == 0)
              continue;
            PCU_ALWAYS_ASSERT(BL_growthRegionsAndLayerFaces
                      ((pRegion)seedRegion,growthRegion,growthLayerFace,Layer_Entity) == 1);
            for (int iglf = 0; iglf < PList_size(growthLayerFace); iglf++) {
              pFace blFace = (pFace)PList_item(growthLayerFace, iglf);
              double fq = ma::measureElementQuality(m, sf, reinterpret_cast<apf::MeshEntity*> (blFace)); // squared mean ratio
              fq = (fq > 0) ? sqrt(fq) : -sqrt(-fq); // mean ratio
              if (fq < minfq)
                minfq = fq;
            }
            PList_clear(growthRegion);
            PList_clear(growthLayerFace);
          }
        }
      }
      FIter_delete(fIter);
    }
    GFIter_delete(gFIter);
    PList_delete(growthRegion);
    PList_delete(growthLayerFace);
  }
  else
#endif
  {
    lion_oprint(1,"PUMI-based mesh doesn't have BL info. minfq is set to be 1.0\n");
  }

// restore mesh
  counter = 0;
  vit = m->begin(0);
  while ((v = m->iterate(vit))) {
// use the original coordinates
    p[0] = op[counter][0];
    p[1] = op[counter][1];
    p[2] = op[counter][2];
    m->setPoint(v, 0, p);
    counter++;
  }
  m->end(vit);
  PCU_ALWAYS_ASSERT(counter==numnp);

// free op
  for (int i = 0; i < numnp; i++)
    free(op[i]);
  free(op);
}

#ifdef __cplusplus
}
#endif

