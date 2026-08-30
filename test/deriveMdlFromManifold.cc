/*
 * This test builds a minimal (6-tetrahedron) mesh of a unit cube with no
 * geometric model, using apf::construct against a null model, and then
 * checks that apf::deriveMdlFromManifold correctly classifies the mesh
 * boundary onto the 6 model faces described by the input bFaces array.
 *
 * This exercises getFaceIdInRegion/getEdgeIdInFace, which look up global
 * vertex ids stored in the internal "_vert_id" mesh tag. That tag is
 * created with createLongTag/setLongTag, so it must be read back with
 * getLongTag; reading it with getIntTag corrupts the comparison on
 * platforms where sizeof(long) != sizeof(int).
 */
#include <apf.h>
#include <apfConvert.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  pcu::Init(&argc, &argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  gmi_register_null();

  // 6-tet fan decomposition of a unit cube around the (0,6) diagonal,
  // vertices ordered as in mds/apfBox.cc's buildTets.
  apf::Gid conn[6][4] = {
    {0,1,2,6},
    {0,2,3,6},
    {0,3,7,6},
    {0,7,4,6},
    {0,4,5,6},
    {0,5,1,6},
  };
  int nelem = 6;
  int nverts = 8;
  double coords[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
    {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
  };

  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false, &PCUObj);
  apf::GlobalToVert globalToVert;
  apf::NewElements elems =
    apf::construct(m, &conn[0][0], nelem, apf::Mesh::TET, globalToVert);
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, &coords[0][0], nverts, globalToVert);
  m->verify();

  std::map<int, apf::MeshEntity*> globalToRegion;
  for (int i = 0; i < nelem; ++i)
    globalToRegion[i] = elems[i];

  bool isModelVert[8];
  for (int i = 0; i < nverts; ++i)
    isModelVert[i] = true;

  // [model_face_tag, adj_region_tag, global_vtx_id_1, global_vtx_id_2,
  //  global_vtx_id_3], two triangles per cube face.
  int nBFaces = 12;
  int bFaces[12][5] = {
    {1, 0, 0, 1, 2},
    {1, 1, 0, 2, 3},
    {2, 0, 1, 2, 6},
    {2, 5, 1, 5, 6},
    {3, 2, 0, 3, 7},
    {3, 3, 0, 4, 7},
    {4, 5, 0, 1, 5},
    {4, 4, 0, 4, 5},
    {5, 1, 2, 3, 6},
    {5, 2, 3, 6, 7},
    {6, 4, 4, 5, 6},
    {6, 3, 4, 6, 7},
  };

  apf::deriveMdlFromManifold(m, isModelVert, nBFaces, bFaces,
                              globalToVert, globalToRegion);

  // Every boundary face must be classified onto the model face tag given
  // in bFaces, and every region must remain on the default model region.
  for (int i = 0; i < nBFaces; ++i) {
    apf::MeshEntity* region = globalToRegion[bFaces[i][1]];
    apf::Downward faces;
    m->getDownward(region, 2, faces);
    apf::MeshEntity* face = NULL;
    apf::Downward verts;
    for (int j = 0; j < 4 && !face; ++j) {
      m->getDownward(faces[j], 0, verts);
      int matched = 0;
      for (int k = 0; k < 3; ++k)
        for (int l = 2; l < 5; ++l)
          if (verts[k] == globalToVert[bFaces[i][l]])
            ++matched;
      if (matched == 3)
        face = faces[j];
    }
    PCU_ALWAYS_ASSERT(face);
    apf::ModelEntity* me = m->toModel(face);
    PCU_ALWAYS_ASSERT(m->getModelType(me) == 2);
    PCU_ALWAYS_ASSERT(m->getModelTag(me) == bFaces[i][0]);
  }

  for (int i = 0; i < nelem; ++i) {
    apf::ModelEntity* me = m->toModel(globalToRegion[i]);
    PCU_ALWAYS_ASSERT(m->getModelType(me) == 3);
    PCU_ALWAYS_ASSERT(m->getModelTag(me) == 0);
  }

  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
  return 0;
}
