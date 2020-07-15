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
#include "apfPolyBasis1D.h"
#include <mth.h>
#include <mth_def.h>
#include <mthQR.h>
#include <pcu_util.h>
#include <PCU.h>

#include <iostream>
using namespace std;

namespace apf {

static unsigned const MAX_ND_ORDER = 10;

// For L2Shapes there are only internal nodes
static inline int countTriNodes(int P)
{
  return (P+1)*(P+2)/2;
}
// For L2Shapes there are only internal nodes
static inline int countTetNodes(int P)
{
  return (P+1)*(P+2)*(P+3)/6;
}

static void computeTriangleTi(
    int P, /*order*/
    mth::Matrix<double>& Q, /*Q in QR factorization of Ti*/
    mth::Matrix<double>& R) /*R in QR factorization of Ti*/
{
  int non = countTriNodes(P);

  apf::NewArray<double> op;
  getOpenPoints(P, op);


  const int p = P;
  apf::NewArray<double> shape_x(p+1);
  apf::NewArray<double> shape_y(p+1);
  apf::NewArray<double> shape_l(p+1);

  apf::DynamicArray<apf::Vector3> nodes (non);

  int o = 0;
  for (int j = 0; j <= p; j++) {
    for (int i = 0; i + j <= p; i++)
      {
        double w = op[i] + op[j] + op[p-i-j];
        nodes[o][0] = op[i]/w;  nodes[o][1] = op[j]/w;  nodes[o][2] = 0.;
        o++;
      }
  }

  // Populate T
  mth::Matrix<double> T(non,non);
  for (int m = 0; m < non; m++)
  {
    o = 0;

    double x = nodes[m][0]; double y = nodes[m][1];

    getChebyshevT(p, x, &shape_x[0]);
    getChebyshevT(p, y, &shape_y[0]);
    getChebyshevT(p, 1. - x - y, &shape_l[0]);

    for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      	T(o++, m) = shape_x[i]*shape_y[j]*shape_l[p-i-j];
  }
  mth::decomposeQR(T, Q, R);
}

static void computeTetTi(
    int P, /*order*/
    mth::Matrix<double>& Q, /*Q in QR factorization of Ti*/
    mth::Matrix<double>& R) /*R in QR factorization of Ti*/
{
  int non = countTetNodes(P);

  apf::NewArray<double> op;
  getOpenPoints(P, op);

  const int p = P;

  apf::NewArray<double> shape_x(p+1);
  apf::NewArray<double> shape_y(p+1);
  apf::NewArray<double> shape_z(p+1);
  apf::NewArray<double> shape_l(p+1);

  apf::DynamicArray<apf::Vector3> nodes (non);
  nodes.setSize(non);

  int o = 0;
  // Region loops to get nodes and dof2tk for regions
  for (int k = 0; k <= p; k++) {
    for (int j = 0; j + k <= p; j++) {
      for (int i = 0; i + j + k <= p; i++) {
        double w = op[i] + op[j] + op[k] + op[p-i-j-k];
        nodes[o][0] = op[i]/w;  nodes[o][1] = op[j]/w;  nodes[o][2] = op[k]/w;
        o++;
      }
    }
  }

  // Populate T
  mth::Matrix<double> T(non,non);
  for (int m = 0; m < non; m++)
  {
    o = 0;
    double x = nodes[m][0]; double y = nodes[m][1]; double z = nodes[m][2];

    getChebyshevT(p, x, &shape_x[0]);
    getChebyshevT(p, y, &shape_y[0]);
    getChebyshevT(p, z, &shape_z[0]);
    getChebyshevT(p, 1. - x - y - z, &shape_l[0]);

    for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
        for (int i = 0; i + j + k <= p; i++)
          T(o++, m) = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[p-i-j-k];
	}
  mth::decomposeQR(T, Q, R);
}

static void getTi(
    int P,
    int type,
    mth::Matrix<double>& Q,
    mth::Matrix<double>& R)
{

  bool cond = (type == apf::Mesh::TRIANGLE || type == apf::Mesh::TET);
  PCU_ALWAYS_ASSERT_VERBOSE(cond,
      "type should be either apf::Mesh::TRIANGLE or apf::Mesh::TET!");

  static apf::NewArray<double> transformQ[apf::Mesh::TYPES][MAX_ND_ORDER+1];
  static apf::NewArray<double> transformR[apf::Mesh::TYPES][MAX_ND_ORDER+1];
  int n = type == apf::Mesh::TRIANGLE ? countTriNodes(P) : countTetNodes(P);

  // get the transform matrices if the are not already computed
  if (!transformQ[type][P].allocated()) {
    mth::Matrix<double> LQ(n,n);
    mth::Matrix<double> LR(n,n);
    type == apf::Mesh::TRIANGLE ?
    	  computeTriangleTi(P, LQ, LR) : computeTetTi(P, LQ, LR);

    transformQ[type][P].allocate(n*n);
    transformR[type][P].allocate(n*n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	      transformQ[type][P][i*n+j] = LQ(i,j);
      	transformR[type][P][i*n+j] = LR(i,j);
      }
    }
  }

  // set Q and R using transformQ and transformR
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Q(i,j) = transformQ[type][P][i*n+j];
      R(i,j) = transformR[type][P][i*n+j];
    }
  }
}

template<int P>
class L2ShapeTri: public FieldShape {
  public:
    L2ShapeTri()
    {
      std::stringstream ss;
      ss << "L2ShapeTri" << P;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const { return name.c_str(); }
    bool isVectorShape() {return false;}
    class Triangle : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<double>& shapes) const
      {
        const int p = P;

      	apf::NewArray<double> shape_x(p+1);
      	apf::NewArray<double> shape_y(p+1);
      	apf::NewArray<double> shape_l(p+1);

        int dof = countNodes();
        mth::Matrix<double> u(dof, 1);

        double x = xi[0]; double y = xi[1];

        getChebyshevT(p, x, &shape_x[0]);
        getChebyshevT(p, y, &shape_y[0]);
        getChebyshevT(p, 1. - x - y, &shape_l[0]);

        int n = 0;
        for (int j = 0; j <= p; j++)
          for (int i = 0; i + j <= p; i++)
            u(n++, 0) = shape_x[i]*shape_y[j]*shape_l[p-i-j];

	      mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TRIANGLE, Q, R);

        mth::Matrix<double> S(dof, 1);
      	for(int i = 0; i < 1; i++) // S = Ti * u
        {
          mth::Vector<double> B (dof);
          mth::Vector<double> X (dof);
          for (int j = 0; j < dof; j++) B[j] = u(j,i); // populate b in QR x = b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(j,i) = X[j]; // populate S with x
        }

        shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  shapes[i] = S(i,0);
        }
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for L2ShapeTri. Aborting()!");
      }
      int countNodes() const {return countTriNodes(P);}
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getVectorValues not implemented \
      	    for L2ShapeTri. Try getValues. Aborting()!");
      }
    };
    EntityShape* getEntityShape(int type)
    {
    	PCU_ALWAYS_ASSERT_VERBOSE(type == Mesh::TRIANGLE,
    			"L2ShapeTri only has entity shapes for TRIANGLEs");
      static Triangle tri;
      return &tri;
    }
    // For the following to member functions we only need to
    // consider the interior nodes, i.e.,
    // Faces: no need to count the nodes associated with bounding edges
    // Tets: no need to count the nodes associated with bounding edges/faces
    bool hasNodesIn(int dimension)
    {
      if (dimension == Mesh::typeDimension[Mesh::TRIANGLE])
      	return true;
      return false;
    }
    int countNodesOn(int type)
    {
      if (type == apf::Mesh::TRIANGLE) return countTriNodes(P);
      return 0;
    }
    int getOrder() {return P;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
      PCU_ALWAYS_ASSERT_VERBOSE(type == Mesh::TRIANGLE,
      	  "getNodeXi for L2ShapeTri can be called only for TRIANGLEs");
      apf::NewArray<double> op;
      getOpenPoints(P, op);
      int c = 0;
      for (int j = 0; j <= P; j++) {
      	for (int i = 0; i + j <= P; i++) {
      	  if (node == c) {
      	    double w = op[i] + op[j] + op[P-i-j];
      	    xi = Vector3( op[i]/w, op[j]/w, 0. );
      	    return;
      	  }
      	  else
      	    c++;
      	}
      }
    }
  private:
    std::string name;
};

template<int P>
class L2ShapeTet: public FieldShape {
  public:
    L2ShapeTet()
    {
      std::stringstream ss;
      ss << "L2ShapeTet_" << P;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const { return name.c_str(); }
    bool isVectorShape() {return false;}
    class Tetrahedron : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<double>& shapes) const
      {
        const int p = P;

      	apf::NewArray<double> shape_x(p+1);
      	apf::NewArray<double> shape_y(p+1);
      	apf::NewArray<double> shape_z(p+1);
      	apf::NewArray<double> shape_l(p+1);

        int dof = countNodes();
        mth::Matrix<double> u(dof, 1);

        double x = xi[0]; double y = xi[1]; double z = xi[2];

        getChebyshevT(p, x, &shape_x[0]);
        getChebyshevT(p, y, &shape_y[0]);
        getChebyshevT(p, z, &shape_z[0]);
        getChebyshevT(p, 1. - x - y - z, &shape_l[0]);

        int n = 0;
        for (int k = 0; k <= p; k++)
          for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            	u(n++, 0) = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[p-i-j-k];


      	mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TET, Q, R);

        mth::Matrix<double> S(dof, 1);
      	for(int i = 0; i < 1; i++) // S = Ti * u
        {
          mth::Vector<double> B (dof);
          mth::Vector<double> X (dof);
          for (int j = 0; j < dof; j++) B[j] = u(j,i); // populate b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(j,i) = X[j]; // populate S with x
        }

        shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  shapes[i] = S(i,0);
        }
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for L2ShapeTet. Aborting()!");
      }
      int countNodes() const {return countTetNodes(P);}
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& /*xi*/, apf::NewArray<apf::Vector3>& /*shapes*/) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getVectorValues not implemented for \
      	    L2ShapeTet. Try getValues. Aborting()!");
      }
    };
    EntityShape* getEntityShape(int type)
    {
      PCU_ALWAYS_ASSERT_VERBOSE(type == Mesh::TET,
      	  "L2ShapeTet only has entity shapes for TETs");
      static Tetrahedron tet;
      return &tet;
    }
    // For the following to member functions we only need to
    // consider the interior nodes, i.e.,
    // Faces: no need to count the nodes associated with bounding edges
    // Tets: no need to count the nodes associated with bounding edges/faces
    bool hasNodesIn(int dimension)
    { 
    	if (dimension == Mesh::typeDimension[Mesh::TET])
    		return true;
      return false;
    }
    int countNodesOn(int type)
    {
      if (type == apf::Mesh::TET) return countTetNodes(P);
      return 0;
    }
    int getOrder() {return P;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
      PCU_ALWAYS_ASSERT_VERBOSE(type == Mesh::TET,
      	  "getNodeXi for L2ShapeTet can be called only for TETs");
      apf::NewArray<double> op;
      getOpenPoints(P, op);
      int c = 0;
      for (int k = 0; k <= P; k++) {
      	for (int j = 0; j + k <= P; j++) {
      	  for (int i = 0; i + j + k <= P; i++) {
      	    if( node == c) {
      	      double w = op[i] + op[j] + op[k] + op[P-i-j-k];
      	      xi = Vector3( op[i]/w, op[j]/w,  op[k]/w );
      	      return;
      	    }
      	    else
      	      c++;
      	  }
      	}
      }
    }
  private:
    std::string name;
};


static apf::FieldShape* getL2ShapeTri(int order)
{
  PCU_ALWAYS_ASSERT_VERBOSE(order >= 0,
      "order is expected to be bigger than or equal to 0!");
  PCU_ALWAYS_ASSERT_VERBOSE(order <= 10,
      "order is expected to be less than or equal to 10!");
  static L2ShapeTri<0>  l2_0;
  static L2ShapeTri<1>  l2_1;
  static L2ShapeTri<2>  l2_2;
  static L2ShapeTri<3>  l2_3;
  static L2ShapeTri<4>  l2_4;
  static L2ShapeTri<5>  l2_5;
  static L2ShapeTri<6>  l2_6;
  static L2ShapeTri<7>  l2_7;
  static L2ShapeTri<8>  l2_8;
  static L2ShapeTri<9>  l2_9;
  static L2ShapeTri<10> l2_10;
  static FieldShape* const l2Shapes[11] = {&l2_0,
                                           &l2_1,
                                           &l2_2,
                                           &l2_3,
                                           &l2_4,
                                           &l2_5,
                                           &l2_6,
                                           &l2_7,
                                           &l2_8,
                                           &l2_9,
                                           &l2_10};
  return l2Shapes[order];
}

static apf::FieldShape* getL2ShapeTet(int order)
{
  PCU_ALWAYS_ASSERT_VERBOSE(order >= 0,
      "order is expected to be bigger than or equal to 0!");
  PCU_ALWAYS_ASSERT_VERBOSE(order <= 10,
      "order is expected to be less than or equal to 10!");
  static L2ShapeTet<0>  l2_0;
  static L2ShapeTet<1>  l2_1;
  static L2ShapeTet<2>  l2_2;
  static L2ShapeTet<3>  l2_3;
  static L2ShapeTet<4>  l2_4;
  static L2ShapeTet<5>  l2_5;
  static L2ShapeTet<6>  l2_6;
  static L2ShapeTet<7>  l2_7;
  static L2ShapeTet<8>  l2_8;
  static L2ShapeTet<9>  l2_9;
  static L2ShapeTet<10> l2_10;
  static FieldShape* const l2Shapes[11] = {&l2_0,
                                           &l2_1,
                                           &l2_2,
                                           &l2_3,
                                           &l2_4,
                                           &l2_5,
                                           &l2_6,
                                           &l2_7,
                                           &l2_8,
                                           &l2_9,
                                           &l2_10};
  return l2Shapes[order];
}


apf::FieldShape* getL2Shape(int order, int type)
{
  if (type == Mesh::TRIANGLE)
    return getL2ShapeTri(order);
  else if (type == Mesh::TET)
    return getL2ShapeTet(order);
  else
    PCU_ALWAYS_ASSERT_VERBOSE(0,
    	"L2Shapes are only implemented for tris and tets");
}

void projectL2Field(Field* to, Field* from)
{
  // checks on the from field
  // checks on the to field
  apf::FieldShape* tShape = getShape(to);
  std::string      tName  = tShape->getName();
  int              tOrder = tShape->getOrder();
  PCU_ALWAYS_ASSERT_VERBOSE((tName == std::string("Linear")) && (tOrder == 1),
                "The to field needs to be 1st order Lagrange!");

  Mesh* m = getMesh(from);
  // auxiliary count fields
  Field* count = createField(m, "counter", SCALAR, getLagrange(1));
  double xis[4][3] = {{0., 0., 0.},
                      {1., 0., 0.},
                      {0., 1., 0.},
                      {0., 0., 1.}};
  // zero out the fields
  zeroField(to);
  zeroField(count);

  int nc = countComponents(to);
  NewArray<double> atXi(nc);
  NewArray<double> currentVal(nc);
  NewArray<double> sum(nc);

  MeshEntity* e;
  MeshIterator* it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    MeshElement* me = createMeshElement(m, e);
    Element* el = createElement(from, me);
    MeshEntity* dvs[4];
    m->getDownward(e, 0, dvs);
    for (int i=0; i<4; i++) {
      getComponents(el, Vector3(xis[i]), &(atXi[0]));
      getComponents(to, dvs[i], 0, &(currentVal[0]));
      for (int j = 0; j < nc; j++) {
        currentVal[j] += atXi[j];
      }
      double currentCount = getScalar(count, dvs[i], 0);
      currentCount += 1.;
      setComponents(to, dvs[i], 0, &(currentVal[0]));
      setScalar(count, dvs[i], 0, currentCount);
    }
    destroyElement(el);
    destroyMeshElement(me);
  }
  m->end(it);

  // take care of entities on part boundary
  accumulate(to);
  accumulate(count);

  it = m->begin(0);
  while( (e = m->iterate(it)) ) {
    getComponents(to, e, 0, &(sum[0]));
    int cnt = getScalar(count, e, 0);
    for (int i = 0; i < nc; i++) {
      sum[i] /= cnt;
    }
    setComponents(to, e, 0, &(sum[0]));
  }
  m->end(it);

  // take care of entities on part boundary
  synchronize(to);

  m->removeField(count);
  destroyField(count);
}

}
