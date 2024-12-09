#include <ma.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <stdlib.h>

#define LEFT_EDGE 0
#define TOP_EDGE 1
#define RIGHT_EDGE 2
#define BOTTOM_EDGE 3
#define TL_VTX 4
#define TR_VTX 5
#define BR_VTX 6
#define BL_VTX 7
#define FACE 8

typedef std::vector<apf::MeshEntity*> EntVec;

void createFace(apf::Mesh2* m, int mdlDim, int mdlId,
    apf::MeshEntity* v0, apf::MeshEntity* v1, apf::MeshEntity* v2) {
  apf::ModelEntity* c = m->findModelEntity(mdlDim, mdlId);
  PCU_ALWAYS_ASSERT(c);
  apf::MeshEntity* faceVerts[3] = {v0,v1,v2};
  apf::buildElement(m, c, apf::Mesh::TRIANGLE, faceVerts);
}

void createEdge(apf::Mesh2* m, int mdlDim, int mdlId,
    apf::MeshEntity* v0, apf::MeshEntity* v1) {
  apf::ModelEntity* c = m->findModelEntity(mdlDim, mdlId);
  PCU_ALWAYS_ASSERT(c);
  apf::MeshEntity* edgeVerts[2] = {v0,v1};
  apf::buildElement(m, c, apf::Mesh::EDGE, edgeVerts);
}

void createVtx(apf::Mesh2* m, int mdlDim, int mdlId, apf::Vector3 pt, EntVec& verts) {
  apf::ModelEntity* c = m->findModelEntity(mdlDim, mdlId);
  PCU_ALWAYS_ASSERT(c);
  apf::MeshEntity* vert = m->createVert(c);
  PCU_ALWAYS_ASSERT(vert);
  m->setPoint(vert, 0, pt);
  verts.push_back(vert);
}

apf::Mesh2* createTriMesh(pcu::PCU *PCUObj) {
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(g, 2, false, PCUObj);

  apf::Vector3 p0( 0.0, 1.0, 0.0);
  apf::Vector3 p1( 0.5, 1.0, 0.0);
  apf::Vector3 p2( 1.0, 1.0, 0.0);
  apf::Vector3 p3( 0.0, 0.5, 0.0);
  apf::Vector3 p4( 0.5, 0.5, 0.0);
  apf::Vector3 p5( 1.0, 0.5, 0.0);
  apf::Vector3 p6( 0.0, 0.0, 0.0);
  apf::Vector3 p7( 0.5, 0.0, 0.0);
  apf::Vector3 p8( 1.0, 0.0, 0.0);

  EntVec verts;
  createVtx(mesh,0,TL_VTX,p0,verts);
  createVtx(mesh,1,TOP_EDGE,p1,verts);
  createVtx(mesh,0,TR_VTX,p2,verts);
  createVtx(mesh,1,LEFT_EDGE,p3,verts);
  createVtx(mesh,2,FACE,p4,verts);
  createVtx(mesh,1,RIGHT_EDGE,p5,verts);
  createVtx(mesh,0,BL_VTX,p6,verts);
  createVtx(mesh,1,BOTTOM_EDGE,p7,verts);
  createVtx(mesh,0,BR_VTX,p8,verts);

  createEdge(mesh,1,TOP_EDGE,verts[0],verts[1]);
  createEdge(mesh,1,TOP_EDGE,verts[1],verts[2]);
  createEdge(mesh,1,LEFT_EDGE,verts[0],verts[3]);
  createEdge(mesh,1,LEFT_EDGE,verts[3],verts[6]);
  createEdge(mesh,1,RIGHT_EDGE,verts[2],verts[5]);
  createEdge(mesh,1,RIGHT_EDGE,verts[5],verts[8]);
  createEdge(mesh,1,BOTTOM_EDGE,verts[6],verts[7]);
  createEdge(mesh,1,BOTTOM_EDGE,verts[7],verts[8]);

  createEdge(mesh,2,FACE,verts[1],verts[4]);
  createEdge(mesh,2,FACE,verts[3],verts[4]);
  createEdge(mesh,2,FACE,verts[4],verts[5]);
  createEdge(mesh,2,FACE,verts[4],verts[7]);

  createEdge(mesh,2,FACE,verts[1],verts[3]);
  createEdge(mesh,2,FACE,verts[2],verts[4]);
  createEdge(mesh,2,FACE,verts[4],verts[6]);
  createEdge(mesh,2,FACE,verts[5],verts[7]);

  createFace(mesh,2,FACE,verts[0],verts[1],verts[3]);
  createFace(mesh,2,FACE,verts[1],verts[4],verts[3]);
  createFace(mesh,2,FACE,verts[1],verts[2],verts[4]);
  createFace(mesh,2,FACE,verts[2],verts[5],verts[4]);
  createFace(mesh,2,FACE,verts[3],verts[4],verts[6]);
  createFace(mesh,2,FACE,verts[4],verts[7],verts[6]);
  createFace(mesh,2,FACE,verts[4],verts[5],verts[7]);
  createFace(mesh,2,FACE,verts[5],verts[8],verts[7]);

  mesh->acceptChanges();
  return mesh;
}

void printClassCounts(apf::Mesh* m) {
  typedef std::map<int,int> mii;
  for(int dim=0; dim<3; dim++) {
    mii meshToMdl;
    apf::MeshIterator* it = m->begin(dim);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::ModelEntity* c = m->toModel(e);
      int mdlTag = m->getModelTag(c);
      meshToMdl[mdlTag]++;
    }
    m->end(it);
    printf("mesh dimension %d\n", dim);
    for(mii::iterator miit = meshToMdl.begin(); 
        miit != meshToMdl.end();
        miit++) {
      printf("mdl id %d has %d mesh ents classified\n", miit->first, miit->second);
    }
  }

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj =  pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);

  gmi_register_null();
  ma::Mesh* m = createTriMesh(&pcu_obj);
  m->verify(); // make sure the mesh is valid!
  gmi_model* g = m->getModel();
  gmi_write_dmg(g,"model.dmg");
  m->writeNative("beforeAdapt.smb");

  printClassCounts(m);

  const ma::Input* in = ma::configureUniformRefine(m, 1);
  ma::adapt(in);

  printClassCounts(m);

  m->writeNative("afterAdapt.smb");

  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}
