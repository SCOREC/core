/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfFieldOf.h"
#include "apfElement.h"
#include "apfVectorElement.h"
#include <cassert>

#include <iostream>

namespace apf {

static const double c0 = -2.449489742783178;
static const double c1 = -3.162277660168379;

static double sign(bool should_flip) {
  return should_flip ? -1.0 : 1.0;
}

class HVertex : public EntityShape {
  public:
    void getValues(
        Mesh*, MeshEntity*, Vector3 const&, NewArray<double>& N) const {
      N.allocate(1);
      N[0] = 1.0;
    }
    void getLocalGradients(
        Mesh*, MeshEntity*, Vector3 const&, NewArray<Vector3>&) const {}
    int countNodes() const { return 1; }
};

class HEdge2 : public EntityShape {
  public:
    void getValues(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<double>& N) const {
      N.allocate(3);
      N[0] = 0.5 * (1.0 - xi[0]);
      N[1] = 0.5 * (1.0 + xi[0]);
      N[2] = 0.25 * c0 * (1.0 - xi[0]*xi[0]);
    }
    void getLocalGradients(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<Vector3>& dN) const {
      dN.allocate(3);
      dN[0] = Vector3(-0.5, 0.0, 0.0);
      dN[1] = Vector3(0.5, 0.0, 0.0);
      dN[2] = Vector3(-0.5 * c0 * xi[0], 0.0, 0.0);
    }
    int countNodes() const { return 3; }
};

class HEdge3 : public EntityShape {
  public:
    void getValues(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<double>& N) const {
      N.allocate(4);
      N[0] = 0.5 * (1.0 - xi[0]);
      N[1] = 0.5 * (1.0 + xi[0]);
      N[2] = 0.25 * c0 * (1.0 - xi[0]*xi[0]);
      N[3] = 0.25 * c0 * (1.0 - xi[0]*xi[0])*xi[0];
    }
    void getLocalGradients(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<Vector3>& dN) const {
      dN.allocate(4);
      dN[0] = Vector3(-0.5, 0.0, 0.0);
      dN[1] = Vector3(0.5, 0.0, 0.0);
      dN[2] = Vector3(-0.5 * c0 * xi[0], 0.0, 0.0);
      dN[3] = Vector3(0.25 * c0 * (1.0 - 3.0*xi[0]*xi[0]), 0.0, 0.0);
    }
    int countNodes() const { return 4; }
};

class HTriangle2 : public EntityShape {
  public:
    void getValues(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<double>& N) const {
      N.allocate(6);

      /* vertices */
      N[0] = 1.0 - xi[0] - xi[1];
      N[1] = xi[0];
      N[2] = xi[1];

      /* edges */
      N[3] = c0 * N[0] * N[1];
      N[4] = c0 * N[1] * N[2];
      N[5] = c0 * N[2] * N[0];
    }
    void getLocalGradients(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<Vector3>& dN) const {
      dN.allocate(6);

      /* vertices */
      dN[0] = Vector3(-1.0, -1.0, 0.0);
      dN[1] = Vector3( 1.0,  0.0, 0.0);
      dN[2] = Vector3( 0.0,  1.0, 0.0);

      /* edges */
      dN[3] = Vector3(c0 * (1.0 - 2.0*xi[0] - xi[1]), -c0 * xi[0], 0.0);
      dN[4] = Vector3(c0 * xi[1], c0 * xi[0], 0.0);
      dN[5] = Vector3(-c0 * xi[1], c0 * (1.0 - xi[0] - 2.0*xi[1]), 0.0);
    }
    int countNodes() const { return 6; }
};

class HTriangle3 : public EntityShape {
  public:
    void getValues(
        Mesh* m, MeshEntity* e, Vector3 const& xi, NewArray<double>& N) const {
      N.allocate(10);

      /* edge orientations */
      int which;
      int rotate;
      bool flip[3];
      double sgn[3];
      apf::MeshEntity* edges[3];
      m->getDownward(e, 1, edges);
      for (int i=0; i < 3; ++i) {
        apf::getAlignment(m, e, edges[i], which, flip[i], rotate);
        sgn[i] = sign(flip[i]);
      }

      /* barycentric coordinates */
      double l0 = 1.0 - xi[0] - xi[1];
      double l1 = xi[0];
      double l2 = xi[1];

      /* vertices */
      N[0] = l0;
      N[1] = l1;
      N[2] = l2;

      /* edge 1 */
      N[3] = c0*l0*l1;
      N[4] = c1*sgn[0]*l0*l1*(l1 - l0);

      /* edge 2 */
      N[5] = c0*l1*l2;
      N[6] = c1*sgn[1]*l1*l2*(l2 - l1);

      /* edge 3 */
      N[7] = c0*l2*l0;
      N[8] = c1*sgn[2]*l2*l0*(l0 - l2);

      /* face */
      N[9] = l0*l1*l2;
    }
    void getLocalGradients(
        Mesh* m, MeshEntity* e, Vector3 const& xi, NewArray<Vector3>& dN) const {
      dN.allocate(10);

      /* edge orientations */
      int which;
      int rotate;
      bool flip[3];
      double sgn[3];
      apf::MeshEntity* edges[3];
      m->getDownward(e, 1, edges);
      for (int i=0; i < 3; ++i) {
        apf::getAlignment(m, e, edges[i], which, flip[i], rotate);
        sgn[i] = sign(flip[i]);
      }

      /* barycentric coordinates */
      double l0 = 1.0 - xi[0] - xi[1];
      double l1 = xi[0];
      double l2 = xi[1];

      /* barycentric gradients */
      apf::Vector3 dl0 = Vector3(-1.0, -1.0, 0.0);
      apf::Vector3 dl1 = Vector3( 1.0,  0.0, 0.0);
      apf::Vector3 dl2 = Vector3( 0.0,  1.0, 0.0);

      /* vertices */
      dN[0] = dl0;
      dN[1] = dl1;
      dN[2] = dl2;

      /* edge 1 */
      dN[3] = (dl0*l1 + dl1*l0)*c0;
      dN[4] = (dl0*l1*(l1-l0) + dl1*l0*(l1-l0) + (dl1-dl0)*l0*l1)*c1*sgn[0];

      /* edge 2 */
      dN[5] = (dl1*l2 + dl2*l1)*c0;
      dN[6] = (dl1*l2*(l2-l1) + dl2*l1*(l2-l1) + (dl2-dl1)*l1*l2)*c1*sgn[1];

      /* edge 3 */
      dN[7] = (dl2*l0 + dl0*l2)*c0;
      dN[8] = (dl2*l0*(l0-l2) + dl0*l2*(l0-l2) + (dl0-dl2)*l2*l0)*c1*sgn[2];

      /* face */
      dN[9] = dl0*l1*l2 + dl1*l0*l2 + dl2*l0*l1;
    }
    void alignSharedNodes(Mesh*, MeshEntity*, MeshEntity*, int order[])
    {
      /* unlike Lagrange shape functions, hierarchic 'modes' do not
         correspond to physical point in space. we want the nodes to
         always come out in terms of dimension-specific polynomial
         order */
      for (int i=0; i < 2; ++i)
        order[i] = i;
    }
    int countNodes() const {return 10;}
};

class HTetrahedron2 : public EntityShape {
  public:
    void getValues(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<double>& N) const {
      N.allocate(10);

      /* vertices */
      N[0] = 1.0 - xi[0] - xi[1] - xi[2];
      N[1] = xi[0];
      N[2] = xi[1];
      N[3] = xi[2];

      /* edges */
      N[4] = c0 * N[0] * N[1];
      N[5] = c0 * N[1] * N[2];
      N[6] = c0 * N[2] * N[0];
      N[7] = c0 * N[0] * N[3];
      N[8] = c0 * N[1] * N[3];
      N[9] = c0 * N[2] * N[3];
    }
    void getLocalGradients(
        Mesh*, MeshEntity*, Vector3 const& xi, NewArray<Vector3>& dN) const {
      dN.allocate(10);

      /* vertices */
      dN[0] = Vector3(-1.0, -1.0, -1.0);
      dN[1] = Vector3( 1.0,  0.0,  0.0);
      dN[2] = Vector3( 0.0,  1.0,  0.0);
      dN[3] = Vector3( 0.0,  0.0,  1.0);

      /* edges */
      dN[4] = Vector3( 1.0-2.0*xi[0]-xi[1]-xi[2], -xi[0], -xi[0] ) * c0;
      dN[5] = Vector3( xi[1], xi[0], 0.0 ) * c0;
      dN[6] = Vector3( -xi[1], 1.0-xi[0]-2.0*xi[1]-xi[2], -xi[1] ) * c0;
      dN[7] = Vector3( -xi[2], -xi[2], 1.0-xi[0]-xi[1]-2.0*xi[2] ) * c0;
      dN[8] = Vector3( xi[2], 0.0, xi[0] ) * c0;
      dN[9] = Vector3( 0.0, xi[2], xi[1] ) * c0;
    }
    int countNodes() const {return 10;}
};

class Hierarchic2 : public FieldShape {
  public:
    const char* getName() const { return "Hierarchic2"; }
    EntityShape* getEntityShape(int type) {
      static HVertex vtx;
      static HEdge2 edge;
      static HTriangle2 tri;
      static HTetrahedron2 tet;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vtx,      //vertex
       &edge,     //edge
       &tri,      //triangle
       NULL,      //quad
       &tet,      //tet
       NULL,      //hex
       NULL,      //prism
       NULL};     //pyramid
      return shapes[type];
    }
    void getNodeXi(int, int, Vector3& xi) {
      xi = Vector3(0,0,0);
    }
    bool hasNodesIn(int dimension) {
      if ( (dimension == 0) || (dimension == 1) )
        return true;
      else
        return false;
    }
    int countNodesOn(int type) {
      if ( (type == Mesh::VERTEX) || (type == Mesh::EDGE) )
        return 1;
      else
        return 0;
    }
    int getOrder() {return 2;}
};

class Hierarchic3 : public FieldShape
{
  public:
    Hierarchic3() {}
    const char* getName() const { return "Hierarchic3"; }
    EntityShape* getEntityShape(int type)
    {
      static HVertex vertex;
      static HEdge3 edge;
      static HTriangle3 tri;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,   //vertex
       &edge,     //edge
       &tri,      //triangle
       NULL,      //quad
       NULL,      //tet
       NULL,      //hex
       NULL,      //prism
       NULL};     //pyramid
      return shapes[type];
    }
    void getNodeXi(int, int, Vector3& xi)
    {
      xi = Vector3(0,0,0);
    }
    bool hasNodesIn(int dimension)
    {
      if ((dimension == 0) ||
          (dimension == 1) ||
          (dimension == 2))
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if ((type == Mesh::VERTEX) || (type == Mesh::TRIANGLE))
        return 1;
      else if (type == Mesh::EDGE)
        return 2;
      else
        return 0;
    }
    int getOrder() {return 3;}
};

FieldShape* getHierarchic(int o)
{
  static Hierarchic2 q;
  static Hierarchic3 c;
  if (o == 1) return getLagrange(o);
  else if (o == 2) return &q;
  else if (o == 3) return &c;
  else return NULL;
}

template <class T>
class Projector : public FieldOp
{
  public:
    Projector(Field* a, Field* b)
    {
      to = static_cast<FieldOf<T>*>(a);
      from = static_cast<FieldOf<T>*>(b);
      mesh = to->getMesh();
      meshElem = 0;
      fromElem = 0;
      int n = to->countComponents();
      data.allocate(n);
      for (int i=0; i < n; ++i)
        data[i] = 0.0;
    }
    bool inEntity(MeshEntity* e)
    {
      meshElem = createMeshElement(mesh, e);
      fromElem = static_cast<ElementOf<T>*>(createElement(from, meshElem));
      return true;
    }
    void atNode(int n)
    {
      int nt = to->countNodesOn(meshElem->getEntity());
      int nf = from->countNodesOn(meshElem->getEntity());
      if ( (nf == 0) || ((nf - nt) < 0))
        setComponents(to, meshElem->getEntity(), n, &data[0]);
      else {
        Vector3 xi;
        to->getShape()->getNodeXi(fromElem->getType(),n,xi);
        T value[1];
        value[0] = fromElem->getValue(xi);
        to->setNodeValue(fromElem->getEntity(),n,value);
      }
    }
    void outEntity()
    {
      destroyMeshElement(meshElem);
      destroyElement(fromElem);
    }
    void run()
    {
      apply(to);
    }
    FieldOf<T>* to;
    FieldOf<T>* from;
    Mesh* mesh;
    VectorElement* meshElem;
    ElementOf<T>* fromElem;
    apf::NewArray<double> data;
};

void projectHierarchicField(Field* to, Field* from)
{
  int ttype = to->getValueType();
  int ftype = from->getValueType();
  assert(ttype == ftype);
  if (ttype == SCALAR) {
    Projector<double> projector(to,from);
    projector.run();
  }
  else if (ttype == VECTOR) {
    Projector<Vector3> projector(to,from);
    projector.run();
  }
  else if (ttype == MATRIX) {
    Projector<Matrix3x3> projector(to,from);
    projector.run();
  }
  else
    fail("projectHierarchicField: unsupported value type");
}

}
