#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>

/* This executable demos how to make a PUMI  mesh using a bottom up approach.
 * The info needed to achieve this:
 * 1) Vert coordinates, and model parametrization
 * 2) Vert classification info (i.e., model dim, model tag)
 * 3) Edge to Vert Connectivity
 * 4) Edge classification info (i.e,. model dim, model tag)
 * 5) Face to Edge Connectivity
 * 6) Face classification info (i.e,. model dim, model tag)
 * 7) Tet to Face Connectivity
 * 8) Tet classification info (i.e,. model dim, model tag)
 *
 * Notes
 * I  ) for this example the above info are provided ad const arrays
 * II ) if model parametrization info is not available
 *      use apf::Vector3(0,0,0), instead
 * III) if model classification info is not available
 *      use 0 is input the corresponding functions
 *  IV) if model info is not available, you might be better of
 *      using the top-down approach instead (refer to construct.cc)
 */

// input data
// each row is x,y,x, p0, p1, p2 of the vertex
const int numV = 9;
const double vert_coords[9][6] = {
{-5.000000e-01, 6.500000e+01, 6.500000e+01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{6.500000e+01, 6.500000e+01, 6.500000e+01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{6.500000e+01, -5.000000e-01, 6.500000e+01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{-5.000000e-01, -5.000000e-01, 6.500000e+01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{-5.000000e-01, 6.500000e+01, -5.000000e-01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{6.500000e+01, 6.500000e+01, -5.000000e-01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{6.500000e+01, -5.000000e-01, -5.000000e-01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{-5.000000e-01, -5.000000e-01, -5.000000e-01, 0.000000e+00, 0.000000e+00, 6.934441e-310},
{2.320018e+01, 5.811620e+01, 3.864902e+01, 0.000000e+00, 0.000000e+00, 6.934441e-310}
};
// each row is dim, tag of the model the vertex is classified on
const int vert_classification[9][2] = {
{0, 22},
{0, 20},
{0, 18},
{0, 16},
{0, 132},
{0, 134},
{0, 136},
{0, 138},
{3, 91}
};
// edge to vertex connectivity
// each row is v0, v1, dim, tag of edge
const int numE = 25;
const int edge_info[25][4] = {
{3, 7, 1, 73},
{2, 6, 1, 72},
{1, 5, 1, 71},
{8, 0, 3, 91},
{0, 4, 1, 70},
{3, 2, 1, 7},
{5, 4, 1, 108},
{8, 4, 3, 91},
{3, 8, 3, 91},
{6, 5, 1, 112},
{1, 0, 1, 13},
{8, 1, 3, 91},
{1, 4, 2, 75},
{0, 2, 2, 1},
{2, 1, 1, 10},
{0, 3, 1, 4},
{7, 6, 1, 116},
{4, 7, 1, 104},
{3, 4, 2, 81},
{3, 6, 2, 79},
{5, 8, 3, 91},
{5, 2, 2, 77},
{3, 5, 3, 91},
{7, 5, 2, 156},
{8, 2, 3, 91}
};
// face to edge connectivity
// each row is e0, e1, e2, dim, tag of face
const int numF = 28;
const int face_info[28][5] = {
{7, 4, 3, 3, 91},
{11, 10, 3, 3, 91},
{12, 11, 7, 3, 91},
{12, 4, 10, 2, 75},
{13, 14, 10, 2, 1},
{8, 3, 15, 3, 91},
{15, 5, 13, 2, 1},
{17, 0, 18, 2, 81},
{18, 8, 7, 3, 91},
{4, 18, 15, 2, 81},
{19, 0, 16, 2, 79},
{19, 1, 5, 2, 79},
{20, 2, 11, 3, 91},
{7, 20, 6, 3, 91},
{12, 2, 6, 2, 75},
{21, 2, 14, 2, 77},
{1, 9, 21, 2, 77},
{20, 8, 22, 3, 91},
{6, 22, 18, 3, 91},
{22, 9, 19, 3, 91},
{23, 9, 16, 2, 156},
{23, 0, 22, 3, 91},
{23, 17, 6, 2, 156},
{24, 20, 21, 3, 91},
{24, 5, 8, 3, 91},
{21, 5, 22, 3, 91},
{14, 11, 24, 3, 91},
{3, 13, 24, 3, 91}
};
// tet to face connectivity
// each row is f0, f1, f2, f3, dim, tag of tet
const int numT = 11;
const int tet_info[11][6] = {
{2, 0, 1, 3, 3, 91},
{8, 9, 5, 0, 3, 91},
{13, 14, 2, 12, 3, 91},
{17, 8, 13, 18, 3, 91},
{20, 10, 19, 21, 3, 91},
{22, 7, 21, 18, 3, 91},
{23, 24, 17, 25, 3, 91},
{15, 12, 23, 26, 3, 91},
{19, 25, 11, 16, 3, 91},
{24, 6, 27, 5, 3, 91},
{26, 4, 1, 27, 3, 91}
};

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);

  PCU_ALWAYS_ASSERT_VERBOSE(PCUObj.Peers() == 1, "Not implemented in parallel!");
  if (argc < 2) {
    printf("USAGE 1 (  no model): %s <outmesh.smb>\n", argv[0]);
    printf("USAGE 2 (with model): %s <outmesh.smb> <model_file\n", argv[0]);
  }
  gmi_register_mesh();
  gmi_register_null();

  const char* outMeshName = argv[1];

  gmi_model* mdl;

  if (argc == 3)
    mdl = gmi_load(argv[2]);
  else
    mdl = gmi_load(".null");

  apf::Mesh2* outMesh = apf::makeEmptyMdsMesh(mdl, 3, false, &PCUObj);

  apf::MeshEntity* verts[numV];
  apf::MeshEntity* edges[numE];
  apf::MeshEntity* faces[numF];

  // construct the verts first
  // Note that each vertex will be added to the verts array in the order
  // it appears in the input above. This way we can easily get the vertex
  // pointer when constructing the edges using edge to vertex connectivity
  for (int i = 0; i < numV; i++) {
    apf::Vector3 coords(vert_coords[i][0],
    	                vert_coords[i][1],
    	                vert_coords[i][2]);
    apf::Vector3 params(vert_coords[i][3],
    	                vert_coords[i][4],
    	                vert_coords[i][5]);
    int modelDim = vert_classification[i][0];
    int modelTag = vert_classification[i][1];

    // find the model entity
    apf::ModelEntity* g = outMesh->findModelEntity(modelDim, modelTag);
    if (g)
      verts[i] = outMesh->createVertex(g, coords, params);
    else
      verts[i] = outMesh->createVertex(0, coords, params);
  }

  // construct the edges next
  // Note that each edge will be added to the edges array in the order
  // it appears in the input above. This way we can easily get the edge
  // pointers when constructing the faces using face to edge connectivity
  for (int i = 0; i < numE; i++) {
    apf::MeshEntity* downVs[2] = {verts[edge_info[i][0]],
                                  verts[edge_info[i][1]]};
    int modelDim = edge_info[i][2];
    int modelTag = edge_info[i][3];

    // find the model entity
    apf::ModelEntity* g = outMesh->findModelEntity(modelDim, modelTag);
    if (g)
      edges[i] = outMesh->createEntity(apf::Mesh::EDGE, g, downVs);
    else
      edges[i] = outMesh->createEntity(apf::Mesh::EDGE, 0, downVs);
  }

  // construct the faces next
  // Note that each face will be added to the faces array in the order
  // it appears in the input above. This way we can easily get the face
  // pointers when constructing the tets using tet to face connectivity
  for (int i = 0; i < numF; i++) {
    apf::MeshEntity* downEs[3] = {edges[face_info[i][0]],
                                  edges[face_info[i][1]],
				  edges[face_info[i][2]]};
    int modelDim = face_info[i][3];
    int modelTag = face_info[i][4];

    // find the model entity
    apf::ModelEntity* g = outMesh->findModelEntity(modelDim, modelTag);
    if (g)
      faces[i] = outMesh->createEntity(apf::Mesh::TRIANGLE, g, downEs);
    else
      faces[i] = outMesh->createEntity(apf::Mesh::TRIANGLE, 0, downEs);
  }

  // finally construct the tets
  for (int i = 0; i < numT; i++) {
    apf::MeshEntity* downFs[4] = {faces[tet_info[i][0]],
                                  faces[tet_info[i][1]],
				  faces[tet_info[i][2]],
				  faces[tet_info[i][3]]};
    int modelDim = tet_info[i][4];
    int modelTag = tet_info[i][5];

    // find the model entity
    apf::ModelEntity* g = outMesh->findModelEntity(modelDim, modelTag);
    if (g)
      outMesh->createEntity(apf::Mesh::TET, g, downFs);
    else
      outMesh->createEntity(apf::Mesh::TET, 0, downFs);
  }

  outMesh->acceptChanges();
  outMesh->verify();
  outMesh->writeNative(outMeshName);


  outMesh->destroyNative();
  apf::destroyMesh(outMesh);
  }
  MPI_Finalize();
}

