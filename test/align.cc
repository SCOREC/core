#include <apfMesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <apf.h>
#include <pcu_util.h>

void testTriEdge(pcu::PCU *PCUObj)
{
  int which, rotate;
  bool flip;
  for(int ed = 0; ed < 6; ++ed){
    gmi_model* model = gmi_load(".null");
    apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false, PCUObj);
    apf::MeshEntity* v[3];
    for (int i = 0; i < 3; ++i)
      v[i] = m->createVert(0);
    apf::MeshEntity* ev[2];
    ev[0] = v[apf::tri_edge_verts[ed % 3][ed >= 3]];
    ev[1] = v[apf::tri_edge_verts[ed % 3][ed <  3]];
    apf::MeshEntity* e =
        apf::buildElement(m, 0, apf::Mesh::EDGE, ev);
    apf::MeshEntity* t =
        apf::buildElement(m, 0, apf::Mesh::TRIANGLE, v);
    apf::getAlignment(m, t, e, which, flip, rotate);
    PCU_ALWAYS_ASSERT(which == (ed % 3));
    PCU_ALWAYS_ASSERT(flip == (ed >= 3));
    PCU_ALWAYS_ASSERT(rotate == 0);
    m->destroyNative();
    apf::destroyMesh(m);
  }
}
void testTetEdge(pcu::PCU *PCUObj)
{
  int which, rotate;
  bool flip;
  for(int ed = 0; ed < 12; ++ed){
    gmi_model* model = gmi_load(".null");
    apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false, PCUObj);
    apf::MeshEntity* v[4];
    for (int i = 0; i < 4; ++i)
      v[i] = m->createVert(0);
    apf::MeshEntity* ev[2];
    ev[0] = v[apf::tet_edge_verts[ed % 6][ed >= 6]];
    ev[1] = v[apf::tet_edge_verts[ed % 6][ed <  6]];
    apf::MeshEntity* e =
        apf::buildElement(m, 0, apf::Mesh::EDGE, ev);
    apf::MeshEntity* t =
        apf::buildElement(m, 0, apf::Mesh::TET, v);
    apf::getAlignment(m, t, e, which, flip, rotate);
    PCU_ALWAYS_ASSERT(which == (ed % 6));
    PCU_ALWAYS_ASSERT(flip == (ed >= 6));
    PCU_ALWAYS_ASSERT(rotate == 0);
    m->destroyNative();
    apf::destroyMesh(m);
  }
}
void testTetTri(pcu::PCU *PCUObj)
{
  int which, rotate;
  bool flip;
  int r[6] = {0,2,1,2,1,0};
  for(int flipped = 0; flipped < 2; ++flipped){
    for(int fa = 0; fa < 12; ++fa){
      gmi_model* model = gmi_load(".null");
      apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false, PCUObj);
      apf::MeshEntity* v[4];
      for (int i = 0; i < 4; ++i)
        v[i] = m->createVert(0);
      apf::MeshEntity* ev[3];
      ev[0] = v[apf::tet_tri_verts[fa % 4][(0+fa/4) % 3]];
      ev[1] = v[apf::tet_tri_verts[fa % 4][(1+fa/4) % 3]];
      ev[2] = v[apf::tet_tri_verts[fa % 4][(2+fa/4) % 3]];
      if(flipped)
        std::swap(ev[1],ev[2]);
      apf::MeshEntity* e =
          apf::buildElement(m, 0, apf::Mesh::TRIANGLE, ev);
      apf::MeshEntity* t =
          apf::buildElement(m, 0, apf::Mesh::TET, v);
      apf::getAlignment(m, t, e, which, flip, rotate);
      PCU_ALWAYS_ASSERT(which == fa % 4);
      PCU_ALWAYS_ASSERT(flip == flipped);
      PCU_ALWAYS_ASSERT(rotate == r[flipped*3+fa/4]);
      m->destroyNative();
      apf::destroyMesh(m);
    }
  }
}
int main()
{
  pcu::Init(0,0);
  {
  pcu::PCU pcu_obj;
  lion_set_verbosity(1);
  gmi_register_null();
  testTriEdge(&pcu_obj);
  testTetEdge(&pcu_obj);
  testTetTri(&pcu_obj);
  }
  pcu::Finalize();
}
