#include <crv.h>
#include <crvBezier.h>
#include <crvSnap.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>
#include <ma.h>

/*
 * This analytic function is an edge in 3D
 */
void vertFunction(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}
// edges go counter clockwise
void edgeFunction(double const p[2], double x[3], void*)
{
  x[0] = p[0]*p[0]*p[0]*p[0];
  x[1] = p[0]*p[0];
  x[2] = p[0];
}
void reparam_zero(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 0;
  to[1] = 0;
}
void reparam_one(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 1;
  to[1] = 0;
}

agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

gmi_model* makeEdgeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edPer = 0;
  double edRan[2] = {0, 1};
  gmi_add_analytic(model, 0, 0, vertFunction, NULL,NULL,NULL);
  gmi_add_analytic(model, 0, 1, vertFunction, NULL,NULL,NULL);
  gmi_ent* ed = gmi_add_analytic(model, 1, 0, edgeFunction, &edPer, &edRan, 0);
  agm_bdry b = add_bdry(model, ed);
  agm_use u0 = add_adj(model, b, 0);
  gmi_add_analytic_reparam(model, u0, reparam_zero, 0);
  agm_use u1 = add_adj(model, b, 1);
  gmi_add_analytic_reparam(model, u1, reparam_one, 0);
  return model;
}
/*
 * Create a mesh with a single edge,
 * Create two edges as an even subdivision,
 * and keep all three around, using the original one
 * to compare correctness of the split.
 *
 */

void testEdgeSubdivision()
{
  for (int o = 1; o <= 6; ++o){

    gmi_model* model = makeEdgeModel();
    apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 1, true);

    apf::ModelEntity* edgeModel = m->findModelEntity(1,0);

    apf::Vector3 points[2] = {apf::Vector3(0,0,0),apf::Vector3(1,1,1)};
    apf::MeshEntity* v[2];
    for (int i = 0; i < 2; ++i)
      v[i] = m->createVertex(m->findModelEntity(0,i),points[i],points[i]);

    apf::MeshEntity* edge = m->createEntity(apf::Mesh::EDGE,edgeModel, v);

    m->acceptChanges();
    m->verify();

    // curve the mesh
    crv::BezierCurver bc(m,o,0,3);
    bc.run();

    apf::Element* elem = apf::createElement(m->getCoordinateField(),edge);
    apf::NewArray<apf::Vector3> nodes;
    apf::NewArray<apf::Vector3> subNodes[2];
    subNodes[0].allocate(o+1);
    subNodes[1].allocate(o+1);
    apf::getVectorNodes(elem,nodes);
    // subdivide the edge's nodes
    crv::subdivideBezierEdge(o,1./4,nodes,subNodes);


    // create the two new edges
    apf::Vector3 p;
    crv::transferParametricOnEdgeSplit(m,edge,1./4,p);

    apf::MeshEntity* v2 = m->createVertex(edgeModel,subNodes[0][o],p);
    apf::MeshEntity* vE[2][2] = {{v[0],v2},{v2,v[1]}};
    apf::MeshEntity* e[2];

    for (int i = 0; i < 2; ++i){
      e[i] = m->createEntity(apf::Mesh::EDGE,edgeModel,vE[i]);
      for (int j = 1; j < o; ++j)
        m->setPoint(e[i],j-1,subNodes[i][j]);
    }
    // compare the two curves to the original one
    apf::Element* elem0 = apf::createElement(m->getCoordinateField(),e[0]);
    apf::Element* elem1 = apf::createElement(m->getCoordinateField(),e[1]);

    apf::Vector3 pt1, pt2, p1, p2;
    for (int i = 0; i <= 100; ++i){
      p1[0] = 0.02*i-1.;
      p2[0] = 0.005*i-1.;
      apf::getVector(elem0,p1,pt1);
      apf::getVector(elem,p2,pt2);

      assert(std::abs((pt2-pt1).getLength() < 1e-15));
      p2[0] = 0.015*i-0.5;
      apf::getVector(elem1,p1,pt1);
      apf::getVector(elem,p2,pt2);
      assert(std::abs((pt2-pt1).getLength() < 1e-15));
    }
    apf::destroyElement(elem);
    apf::destroyElement(elem0);
    apf::destroyElement(elem1);

    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  testEdgeSubdivision();
  PCU_Comm_Free();
  MPI_Finalize();
}

