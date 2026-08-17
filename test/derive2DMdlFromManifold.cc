/*
 * This test builds a minimal (2-triangle) mesh of a unit square with no
 * geometric model, using apf::construct against a null model, and then
 * checks that apf::derive2DMdlFromManifold correctly classifies the mesh
 * boundary onto the 4 model edges described by the input bEdges array.
 *
 * This exercises getEdgeIdInFace, which looks up global vertex ids stored
 * in the internal "_vert_id" mesh tag. That tag is created with
 * createLongTag/setLongTag, so it must be read back with getLongTag;
 * reading it with getIntTag corrupts the comparison on platforms where
 * sizeof(long) != sizeof(int).
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

  // 2-triangle decomposition of a unit square along the (0,2) diagonal.
  apf::Gid conn[2][3] = {
    {0,1,2},
    {0,2,3},
  };
  int nelem = 2;
  int nverts = 4;
  double coords[4][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
  };

  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false, &PCUObj);
  apf::GlobalToVert globalToVert;
  apf::NewElements elems =
    apf::construct(m, &conn[0][0], nelem, apf::Mesh::TRIANGLE, globalToVert);
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, &coords[0][0], nverts, globalToVert);
  m->verify();

  std::map<int, apf::MeshEntity*> globalToFace;
  for (int i = 0; i < nelem; ++i)
    globalToFace[i] = elems[i];

  bool isModelVert[4];
  for (int i = 0; i < nverts; ++i)
    isModelVert[i] = true;

  // [model_edge_tag, adj_face_tag, global_vtx_id_1, global_vtx_id_2],
  // one edge per square side.
  int nBEdges = 4;
  int bEdges[4][4] = {
    {1, 0, 0, 1},
    {2, 0, 1, 2},
    {3, 1, 2, 3},
    {4, 1, 0, 3},
  };

  apf::derive2DMdlFromManifold(m, isModelVert, nBEdges, bEdges,
                                globalToVert, globalToFace);

  // Every boundary edge must be classified onto the model edge tag given
  // in bEdges, and every face must remain on the default model region.
  for (int i = 0; i < nBEdges; ++i) {
    apf::MeshEntity* face = globalToFace[bEdges[i][1]];
    apf::Downward edges;
    m->getDownward(face, 1, edges);
    apf::MeshEntity* edge = NULL;
    apf::Downward verts;
    for (int j = 0; j < 3 && !edge; ++j) {
      m->getDownward(edges[j], 0, verts);
      int matched = 0;
      for (int k = 0; k < 2; ++k)
        for (int l = 2; l < 4; ++l)
          if (verts[k] == globalToVert[bEdges[i][l]])
            ++matched;
      if (matched == 2)
        edge = edges[j];
    }
    PCU_ALWAYS_ASSERT(edge);
    apf::ModelEntity* me = m->toModel(edge);
    PCU_ALWAYS_ASSERT(m->getModelType(me) == 1);
    PCU_ALWAYS_ASSERT(m->getModelTag(me) == bEdges[i][0]);
  }

  for (int i = 0; i < nelem; ++i) {
    apf::ModelEntity* me = m->toModel(globalToFace[i]);
    PCU_ALWAYS_ASSERT(m->getModelType(me) == 2);
    PCU_ALWAYS_ASSERT(m->getModelTag(me) == 0);
  }

  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
  return 0;
}
