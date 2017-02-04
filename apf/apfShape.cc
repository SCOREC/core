/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfIntegrate.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include <cassert>

namespace apf {

static std::map<std::string, FieldShape*> registry;

EntityShape::~EntityShape()
{
}

void EntityShape::alignSharedNodes(Mesh* m,
    MeshEntity* elem, MeshEntity* shared, int order[])
{
  (void)m;
  (void)elem;
  (void)shared;
  (void)order;
  fail("unimplemented alignSharedNodes\n");
}

FieldShape::~FieldShape()
{
}

void FieldShape::getNodeXi(int, int, Vector3&)
{
  fail("unimplemented getNodeXi called");
}

void FieldShape::registerSelf(const char* name_)
{
  std::string name = name_;
  assert(registry.count(name) == 0);
  registry[name] = this;
}

FieldShape* getShapeByName(const char* name)
{
  /* Static variables in functions (which is what
     all the FieldShape objects are) are constructed
     *on the first call to the function*, and they register
     themselves in their constructors.
     If we had an apf::initialize() function, that would be a good
     place to do this, but we don't so we'll do it here.
     Users who have their own FieldShapes should make sure they are
     constructed before this gets called as well. */
  getLagrange(1);
  getSerendipity();
  getConstant(0);
  getIPShape(2,1);
  getVoronoiShape(2,1);
  getIPFitShape(2,1);
  std::string s(name);
  if (registry.count(s))
    return registry[s];
  return 0;
}

class Linear : public FieldShape
{
  public:
    Linear() { registerSelf(apf::Linear::getName()); }
    const char* getName() const { return "Linear"; }
    class Vertex : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1.0;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
        }
        int countNodes() const {return 1;}
    };
    class Edge : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(2);
          values[0] = (1.0-xi[0])/2.0;
          values[1] = (1.0+xi[0])/2.0;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
        }
        int countNodes() const {return 2;}
    };
    class Triangle : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(3);
          values[0] = 1-xi[0]-xi[1];
          values[1] = xi[0];
          values[2] = xi[1];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          grads.allocate(3);
          grads[0] = Vector3(-1,-1,0);
          grads[1] = Vector3( 1, 0,0);
          grads[2] = Vector3( 0, 1,0);
        }
        int countNodes() const {return 3;}
    };
    class Quad : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(4);
          double l0x = (1-xi[0]);
          double l1x = (1+xi[0]);
          double l0y = (1-xi[1]);
          double l1y = (1+xi[1]);
          values[0] = l0x*l0y/4;
          values[1] = l1x*l0y/4;
          values[2] = l1x*l1y/4;
          values[3] = l0x*l1y/4;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          grads.allocate(4);
          double l0x = (1-xi[0]);
          double l1x = (1+xi[0]);
          double l0y = (1-xi[1]);
          double l1y = (1+xi[1]);
          grads[0] = Vector3(-l0y,-l0x,0)/4;
          grads[1] = Vector3( l0y,-l1x,0)/4;
          grads[2] = Vector3( l1y, l1x,0)/4;
          grads[3] = Vector3(-l1y, l0x,0)/4;
        }
        int countNodes() const {return 4;}
    };
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(4);
          values[0] = 1-xi[0]-xi[1]-xi[2];
          values[1] = xi[0];
          values[2] = xi[1];
          values[3] = xi[2];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          grads.allocate(4);
          grads[0] = Vector3(-1,-1,-1);
          grads[1] = Vector3( 1, 0, 0);
          grads[2] = Vector3( 0, 1, 0);
          grads[3] = Vector3( 0, 0, 1);
        }
        int countNodes() const {return 4;}
    };
    class Prism : public EntityShape
    {
      public: //tensor product of triangle and edge
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(6);
          double nt[3];
          nt[0] = 1-xi[0]-xi[1];
          nt[1] = xi[0];
          nt[2] = xi[1];
          double down = (1 - xi[2]) / 2.0;
          double up   = (1 + xi[2]) / 2.0;
          values[0] = nt[0] * down;
          values[1] = nt[1] * down;
          values[2] = nt[2] * down;
          values[3] = nt[0] * up  ;
          values[4] = nt[1] * up  ;
          values[5] = nt[2] * up  ;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          grads.allocate(6);
          double nt[3];
          nt[0] = 1-xi[0]-xi[1];
          nt[1] = xi[0];
          nt[2] = xi[1];
          double const down = (1 - xi[2]) / 2.0;
          double const up   = (1 + xi[2]) / 2.0;
          static Vector3 const tg[3] =
          {Vector3(-1,-1,0),
           Vector3( 1, 0,0),
           Vector3( 0, 1,0)};
          static Vector3 const eg[2] =
          {Vector3(0,0,-0.5),
           Vector3(0,0, 0.5)};
          grads[0] = tg[0] * down + eg[0] * nt[0];
          grads[1] = tg[1] * down + eg[0] * nt[1];
          grads[2] = tg[2] * down + eg[0] * nt[2];
          grads[3] = tg[0] * up   + eg[1] * nt[0];
          grads[4] = tg[1] * up   + eg[1] * nt[1];
          grads[5] = tg[2] * up   + eg[1] * nt[2];
        }
        int countNodes() const {return 6;}
    };
    class Pyramid : public EntityShape
    {
      public: /* degenerate hexahedron */
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(5);
          double l0x = (1 - xi[0]);
          double l1x = (1 + xi[0]);
          double l0y = (1 - xi[1]);
          double l1y = (1 + xi[1]);
          double l0z = (1 - xi[2]);
          double l1z = (1 + xi[2]);
          values[0] = l0x * l0y * l0z / 8;
          values[1] = l1x * l0y * l0z / 8;
          values[2] = l1x * l1y * l0z / 8;
          values[3] = l0x * l1y * l0z / 8;
          values[4] = l1z / 2;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          double l0x = (1 - xi[0]);
          double l1x = (1 + xi[0]);
          double l0y = (1 - xi[1]);
          double l1y = (1 + xi[1]);
          double l0z = (1 - xi[2]);
          grads.allocate(5);
          grads[0] = Vector3(-l0y * l0z, -l0x * l0z, -l0x * l0y) / 8;
          grads[1] = Vector3( l0y * l0z, -l1x * l0z, -l1x * l0y) / 8;
          grads[2] = Vector3( l1y * l0z,  l1x * l0z, -l1x * l1y) / 8;
          grads[3] = Vector3(-l1y * l0z,  l0x * l0z, -l0x * l1y) / 8;
          grads[4] = Vector3(0,0,0.5);
        }
        int countNodes() const {return 5;}
    };
    class Hexahedron : public EntityShape
    {
      public: /* degenerate hexahedron */
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(8);
          double l0x = (1 - xi[0]);
          double l1x = (1 + xi[0]);
          double l0y = (1 - xi[1]);
          double l1y = (1 + xi[1]);
          double l0z = (1 - xi[2]);
          double l1z = (1 + xi[2]);
          values[0] = l0x * l0y * l0z / 8;
          values[1] = l1x * l0y * l0z / 8;
          values[2] = l1x * l1y * l0z / 8;
          values[3] = l0x * l1y * l0z / 8;
          values[4] = l0x * l0y * l1z / 8;
          values[5] = l1x * l0y * l1z / 8;
          values[6] = l1x * l1y * l1z / 8;
          values[7] = l0x * l1y * l1z / 8;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          double l0x = (1 - xi[0]);
          double l1x = (1 + xi[0]);
          double l0y = (1 - xi[1]);
          double l1y = (1 + xi[1]);
          double l0z = (1 - xi[2]);
          double l1z = (1 + xi[2]);
          grads.allocate(8);
          grads[0] = Vector3(-l0y * l0z, -l0x * l0z, -l0x * l0y) / 8;
          grads[1] = Vector3( l0y * l0z, -l1x * l0z, -l1x * l0y) / 8;
          grads[2] = Vector3( l1y * l0z,  l1x * l0z, -l1x * l1y) / 8;
          grads[3] = Vector3(-l1y * l0z,  l0x * l0z, -l0x * l1y) / 8;
          grads[4] = Vector3(-l0y * l1z, -l0x * l1z,  l0x * l0y) / 8;
          grads[5] = Vector3( l0y * l1z, -l1x * l1z,  l1x * l0y) / 8;
          grads[6] = Vector3( l1y * l1z,  l1x * l1z,  l1x * l1y) / 8;
          grads[7] = Vector3(-l1y * l1z,  l0x * l1z,  l0x * l1y) / 8;
        }
        int countNodes() const {return 8;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
      static Edge edge;
      static Triangle triangle;
      static Quad quad;
      static Tetrahedron tet;
      static Prism prism;
      static Pyramid pyramid;
      static Hexahedron hex;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,
       &edge,
       &triangle,
       &quad,
       &tet,
       &hex,
       &prism,
       &pyramid};
      return shapes[type];
    }
    bool hasNodesIn(int dimension)
    {
      if (dimension == 0)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
        return 1;
      else
        return 0;
    }
    int getOrder() {return 1;}
    void getNodeXi(int, int, Vector3& xi)
    {
      xi = Vector3(0,0,0);
    }
};

class QuadraticBase : public FieldShape
{
  public:
    class Edge : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(3);
          values[0] = -xi[0]*(1-xi[0])/2.0;
          values[1] =  xi[0]*(1+xi[0])/2.0;
          values[2] = 1-xi[0]*xi[0];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          grads.allocate(3);
          grads[0] = Vector3((2*xi[0]-1)/2.0,0,0);
          grads[1] = Vector3((2*xi[0]+1)/2.0,0,0);
          grads[2] = Vector3(-2*xi[0],0,0);
        }
        int countNodes() const {return 3;}
    };
    class Triangle : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          double xi2 = 1-xi[0]-xi[1];
          values.allocate(6);
          values[0] = xi2*(2*xi2-1);
          values[1] = xi[0]*(2*xi[0]-1);
          values[2] = xi[1]*(2*xi[1]-1);
          values[3] = 4*xi[0]*xi2;
          values[4] = 4*xi[0]*xi[1];
          values[5] = 4*xi[1]*xi2;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          double xi2 = 1-xi[0]-xi[1];
          grads.allocate(6);
          grads[0] = Vector3(-4*xi2+1,-4*xi2+1,0);
          grads[1] = Vector3(4*xi[0]-1,0,0);
          grads[2] = Vector3(0,4*xi[1]-1,0);
          grads[3] = Vector3(4*(xi2-xi[0]),-4*xi[0],0);
          grads[4] = Vector3(4*xi[1],4*xi[0],0);
          grads[5] = Vector3(-4*xi[1],4*(xi2-xi[1]),0);
        }
        int countNodes() const {return 6;}
    };
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          double xi3 = 1-xi[0]-xi[1]-xi[2];
          values.allocate(10);
          values[0] = xi3*(2*xi3-1);
          values[1] = xi[0]*(2*xi[0]-1);
          values[2] = xi[1]*(2*xi[1]-1);
          values[3] = xi[2]*(2*xi[2]-1);
          values[4] = 4*xi[0]*xi3;
          values[5] = 4*xi[0]*xi[1];
          values[6] = 4*xi[1]*xi3;
          values[7] = 4*xi[2]*xi3;
          values[8] = 4*xi[2]*xi[0];
          values[9] = 4*xi[1]*xi[2];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          double xi3 = 1-xi[0]-xi[1]-xi[2];
          grads.allocate(10);
          double d3 = 1-4*xi3;
          grads[0] = Vector3(d3,d3,d3);
          grads[1] = Vector3(4*xi[0]-1,0,0);
          grads[2] = Vector3(0,4*xi[1]-1,0);
          grads[3] = Vector3(0,0,4*xi[2]-1);
          grads[4] = Vector3(4*xi3-4*xi[0],-4*xi[0],-4*xi[0]);
          grads[5] = Vector3(4*xi[1],4*xi[0],0);
          grads[6] = Vector3(-4*xi[1],4*xi3-4*xi[1],-4*xi[1]);
          grads[7] = Vector3(-4*xi[2],-4*xi[2],4*xi3-4*xi[2]);
          grads[8] = Vector3(4*xi[2],0,4*xi[0]);
          grads[9] = Vector3(0,4*xi[2],4*xi[1]);
        }
        int countNodes() const {return 10;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Linear::Vertex vertex;
      static Edge edge;
      static Triangle triangle;
      static Tetrahedron tet;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,   //vertex
       &edge,     //edge
       &triangle, //triangle
       NULL,     //quad
       &tet,      //tet
       NULL,      //hex
       NULL,      //prism
       NULL};     //pyramid
      return shapes[type];
    }
    int getOrder() {return 2;}
    void getNodeXi(int, int, Vector3& xi)
    {
      /* for vertex nodes, mid-edge nodes,
         and mid-quad nodes, the xi
         coordinate is zero */
      xi = Vector3(0,0,0);
    }
};

class LagrangeQuadratic : public QuadraticBase
{
  public:
    LagrangeQuadratic() { registerSelf(apf::LagrangeQuadratic::getName()); }
    const char* getName() const {return "Lagrange Quadratic";}
    class Quad : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          /*This is ordering of shape functions used below
          *      ^n
          *      |
          *   4--7--3
          *   |  |  |
          *   8  9--6--->e
          *   |     |
          *   1--5--2
          *
          *   indices are just above minus one
          */
          double e = xi[0];
          double n = xi[1];
          values.allocate(9);
          values[0] = (e*n)*(e-1)*(n-1)/4.0;
          values[1] = (e*n)*(e+1)*(n-1)/4.0;
          values[2] = (e*n)*(e+1)*(n+1)/4.0;
          values[3] = (e*n)*(e-1)*(n+1)/4.0;
          values[4] = (1-(e*e))*(n-1)*n/2.0;
          values[5] = (e+1)*(1-(n*n))*e/2.0;
          values[6] = (1-(e*e))*(n+1)*n/2.0;
          values[7] = (e-1)*(1-(n*n))*e/2.0;
          values[8] = (1-(e*e))*(1-(n*n));
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          double e = xi[0];
          double n = xi[1];
          grads.allocate(9);
          grads[0] = Vector3(
              (e - 0.5)*(n*(n-1))/2.0,
              (n - 0.5)*(e*(e-1))/2.0, 0.0);
          grads[1] = Vector3(
              (e + 0.5)*(n*(n-1))/2.0,
              (n - 0.5)*(e*(e+1))/2.0, 0.0);
          grads[2] = Vector3(
              (e + 0.5)*(n*(n+1))/2.0,
              (n + 0.5)*(e*(e+1))/2.0, 0.0);
          grads[3] = Vector3(
              (e - 0.5)*(n*(n+1))/2.0,
              (n + 0.5)*(e*(e-1))/2.0, 0.0);
          grads[4] = Vector3(
              -e*(n*(n-1)),
              (n - 0.5)*(1-(e*e)), 0.0);
          grads[5] = Vector3(
              (e + 0.5)*(1-(n*n)),
              -n*(e*(e+1)), 0.0);
          grads[6] = Vector3(
              -e*(n*(n+1)),
              (n + 0.5)*(1-(e*e)), 0.0);
          grads[7] = Vector3(
              (e - 0.5)*(1-(n*n)),
              -n*(e*(e-1)), 0.0);
          grads[8] = Vector3(
              -2.0*e*(1-(n*n)),
              -2.0*n*(1-(e*e)), 0.0);
        }
        int countNodes() const {return 9;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Quad quad;
      if (type == Mesh::QUAD)
        return &quad;
      return this->QuadraticBase::getEntityShape(type);
    }
    bool hasNodesIn(int dimension)
    {
      if ((dimension == 0)||
          (dimension == 1)||
          (dimension == 2))
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if ((type == Mesh::VERTEX)||
          (type == Mesh::EDGE)||
          (type == Mesh::QUAD))
        return 1;
      else
        return 0;
    }
};

class SerendipityQuadratic : public QuadraticBase
{
  public:
    SerendipityQuadratic() { registerSelf(apf::SerendipityQuadratic::getName()); }
    const char* getName() const {return "Serendipity Quadratic";}
    class Quad : public EntityShape
    {
    public:
      void getValues(Mesh*, MeshEntity*,
          Vector3 const& xi, NewArray<double>& values) const
      {
        values.allocate(8);
        values[0] = (1-xi[0])*(1-xi[1])*(-xi[0] - xi[1] -1)/4.0;
        values[1] = (1+xi[0])*(1-xi[1])*(xi[0] - xi[1] -1)/4.0;
        values[2] = (1+xi[0])*(1+xi[1])*(xi[0] + xi[1] -1)/4.0;
        values[3] = (1-xi[0])*(1+xi[1])*(-xi[0] + xi[1] - 1)/4.0;
        values[4] = (1-xi[0]*xi[0])*(1-xi[1])/2.0;
        values[5] = (1+xi[0])*(1-xi[1]*xi[1])/2.0;
        values[6] = (1-xi[0]*xi[0])*(1+xi[1])/2.0;
        values[7] = (1-xi[0])*(1-xi[1]*xi[1])/2.0;
      }
      void getLocalGradients(Mesh*, MeshEntity*,
          Vector3 const& xi,
          NewArray<Vector3>& grads) const
      {
        grads.allocate(8);
        grads[0] = Vector3(
            xi[0]/2.0 + xi[1]/4.0 - (xi[0]*xi[1])/2.0 - xi[1]*xi[1]/4.0,
            xi[0]/4.0 + xi[1]/2.0 - (xi[0]*xi[1])/2.0 - xi[0]*xi[0]/4.0, 0.0);
        grads[1] = Vector3(
            xi[0]/2.0 - xi[1]/4.0 - (xi[0]*xi[1])/2.0 + xi[1]*xi[1]/4.0,
            xi[1]/2.0 - xi[0]/4.0 + (xi[0]*xi[1])/2.0 - xi[0]*xi[0]/4.0, 0.0);
        grads[2] = Vector3(
            xi[0]/2.0 + xi[1]/4.0 + (xi[0]*xi[1])/2.0 + xi[1]*xi[1]/4.0,
            xi[0]/4.0 + xi[1]/2.0 + (xi[0]*xi[1])/2.0 + xi[0]*xi[0]/4.0, 0.0);
        grads[3] = Vector3(
            xi[0]/2.0 - xi[1]/4.0 + (xi[0]*xi[1])/2.0 - xi[1]*xi[1]/4.0,
            xi[1]/2.0 - xi[0]/4.0 - (xi[0]*xi[1])/2.0 + xi[0]*xi[0]/4.0, 0.0);
        grads[4] = Vector3(xi[0]*xi[1] - xi[0],   xi[0]*xi[0]/2.0 - 0.5, 0.0);
        grads[5] = Vector3(0.5 - xi[1]*xi[1]/2.0,  -xi[1] - xi[0]*xi[1], 0.0);
        grads[6] = Vector3(-xi[0] - xi[0]*xi[1],  0.5 - xi[0]*xi[0]/2.0, 0.0);
        grads[7] = Vector3(xi[1]*xi[1]/2.0 - 0.5,   xi[0]*xi[1] - xi[1], 0.0);
      }
      int countNodes() const {return 8;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Quad quad;
      if (type == Mesh::QUAD)
        return &quad;
      return this->QuadraticBase::getEntityShape(type);
    }
    bool hasNodesIn(int dimension)
    {
      if ((dimension == 0)||
          (dimension == 1))
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if ((type == Mesh::VERTEX)||
          (type == Mesh::EDGE))
        return 1;
      else
        return 0;
    }
};

class LagrangeCubic : public FieldShape
{
  public:
    LagrangeCubic() { registerSelf(apf::LagrangeCubic::getName()); }
    const char* getName() const { return "Lagrange Cubic"; }
    class Vertex : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1.0;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
        }
        int countNodes() const {return 1;}
    };
    class Edge : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& p, NewArray<double>& N) const
        {
          N.allocate(4);
          double xi = p[0];
          N[0] = 9./16.*(1./9.-xi*xi)*(xi-1.);
          N[1] = -9./16.*(1./9.-xi*xi)*(xi+1.);
          N[2] = 27./16.*(1.-xi*xi)*(1./3.-xi);
          N[3] = 27./16.*(1.-xi*xi)*(1./3.+xi);
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& p, NewArray<Vector3>& dN) const
        {
          dN.allocate(4);
          double xi = p[0];
          dN[0] = Vector3(-9./16.*(3.*xi*xi-2.*xi-1./9.), 0, 0);
          dN[1] = Vector3(-9./16.*(-3.*xi*xi-2.*xi+1./9.), 0, 0);
          dN[2] = Vector3(27./16.*(3.*xi*xi-2./3.*xi-1.), 0, 0);
          dN[3] = Vector3(27./16.*(-3.*xi*xi-2./3.*xi+1.), 0, 0);
        }
        int countNodes() const {return 4;}
    };
    class Triangle : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& N) const
        {
          N.allocate(10);
          double l1 = 1.0-xi[0]-xi[1];
          double l2 = xi[0];
          double l3 = xi[1];
          N[0] = 0.5*(3.*l1-1.)*(3.*l1-2.)*l1;
          N[1] = 0.5*(3.*l2-1.)*(3.*l2-2.)*l2;
          N[2] = 0.5*(3.*l3-1.)*(3.*l3-2.)*l3;
          N[3] = 9./2.*l1*l2*(3.*l1-1.);
          N[4] = 9./2.*l1*l2*(3.*l2-1.);
          N[5] = 9./2.*l2*l3*(3.*l2-1.);
          N[6] = 9./2.*l2*l3*(3.*l3-1.);
          N[7] = 9./2.*l3*l1*(3.*l3-1.);
          N[8] = 9./2.*l3*l1*(3.*l1-1.);
          N[9] = 27.*l1*l2*l3;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          dN.allocate(10);
          double l1 = 1.0-xi[0]-xi[1];
          double l2 = xi[0];
          double l3 = xi[1];
          apf::Vector3 gl1(-1,-1,0);
          apf::Vector3 gl2(1,0,0);
          apf::Vector3 gl3(0,1,0);
          dN[0] = gl1*(27./2.*l1*l1-9.*l1+1.);
          dN[1] = gl2*(27./2.*l2*l2-9.*l2+1.);
          dN[2] = gl3*(27./2.*l3*l3-9.*l3+1.);
          dN[3] = (gl2*(3.*l1*l1-l1) + gl1*(6.*l1*l2-l2))*9./2;
          dN[4] = (gl1*(3.*l2*l2-l2) + gl2*(6.*l1*l2-l1))*9./2.;
          dN[5] = (gl3*(3.*l2*l2-l2) + gl2*(6.*l2*l3-l3))*9./2.;
          dN[6] = (gl2*(3.*l3*l3-l3) + gl3*(6.*l2*l3-l2))*9./2.;
          dN[7] = (gl1*(3.*l3*l3-l3) + gl3*(6.*l3*l1-l1))*9./2.;
          dN[8] = (gl3*(3.*l1*l1-l1) + gl1*(6.*l3*l1-l3))*9./2.;
          dN[9] = (gl1*l2*l3+gl2*l1*l3+gl3*l1*l2)*27.;
        }
        int countNodes() const {return 10;}
        void alignSharedNodes(Mesh* m, MeshEntity* elem,
            MeshEntity* edge, int order[]) {
          int which, rotate;
          bool flip;
          getAlignment(m, elem, edge, which, flip, rotate);
          if (! flip) {
            order[0] = 0;
            order[1] = 1;
          }
          else {
            order[0] = 1;
            order[1] = 0;
          }
        }
    };
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& N) const
        {
          N.allocate(20);
          double l0 = 1.0-xi[0]-xi[1]-xi[2];
          double l1 = xi[0];
          double l2 = xi[1];
          double l3 = xi[2];
          N[0] = 0.5*(3.*l0-1.)*(3.*l0-2.)*l0;
          N[1] = 0.5*(3.*l1-1.)*(3.*l1-2.)*l1;
          N[2] = 0.5*(3.*l2-1.)*(3.*l2-2.)*l2;
          N[3] = 0.5*(3.*l3-1.)*(3.*l3-2.)*l3;
          N[4] = 9./2.*l0*l1*(3.*l0-1.);
          N[5] = 9./2.*l0*l1*(3.*l1-1.);
          N[6] = 9./2.*l1*l2*(3.*l1-1.);
          N[7] = 9./2.*l1*l2*(3.*l2-1.);
          N[8] = 9./2.*l2*l0*(3.*l2-1.);
          N[9] = 9./2.*l2*l0*(3.*l0-1.);
          N[10] = 9./2.*l0*l3*(3.*l0-1.);
          N[11] = 9./2.*l0*l3*(3.*l3-1.);
          N[12] = 9./2.*l1*l3*(3.*l1-1.);
          N[13] = 9./2.*l1*l3*(3.*l3-1.);
          N[14] = 9./2.*l2*l3*(3.*l2-1.);
          N[15] = 9./2.*l2*l3*(3.*l3-1.);
          N[16] = 27.*l0*l1*l2;
          N[17] = 27.*l0*l1*l3;
          N[18] = 27.*l1*l2*l3;
          N[19] = 27.*l0*l2*l3;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          dN.allocate(20);
          double l0 = 1.0-xi[0]-xi[1]-xi[2];
          double l1 = xi[0];
          double l2 = xi[1];
          double l3 = xi[2];
          apf::Vector3 gl0(-1,-1,-1);
          apf::Vector3 gl1(1,0,0);
          apf::Vector3 gl2(0,1,0);
          apf::Vector3 gl3(0,0,1);
          dN[0] = gl0*(27./2.*l0*l0-9.*l0+1.);
          dN[1] = gl1*(27./2.*l1*l1-9.*l1+1.);
          dN[2] = gl2*(27./2.*l2*l2-9.*l2+1.);
          dN[3] = gl3*(27./2.*l3*l3-9.*l3+1.);
          dN[4] = (gl1*(3.*l0*l0-l0) + gl0*(6.*l0*l1-l1))*9./2.;
          dN[5] = (gl0*(3.*l1*l1-l1) + gl1*(6.*l0*l1-l0))*9./2.;
          dN[6] = (gl2*(3.*l1*l1-l1) + gl1*(6.*l1*l2-l2))*9./2.;
          dN[7] = (gl1*(3.*l2*l2-l2) + gl2*(6.*l1*l2-l1))*9./2.;
          dN[8] = (gl0*(3.*l2*l2-l2) + gl2*(6.*l2*l0-l0))*9./2.;
          dN[9] = (gl2*(3.*l0*l0-l0) + gl0*(6.*l2*l0-l2))*9./2.;
          dN[10] = (gl3*(3.*l0*l0-l0) + gl0*(6.*l0*l3-l3))*9./2.;
          dN[11] = (gl0*(3.*l3*l3-l3) + gl3*(6.*l0*l3-l0))*9./2.;
          dN[12] = (gl3*(3.*l1*l1-l1) + gl1*(6.*l1*l3-l3))*9./2.;
          dN[13] = (gl1*(3.*l3*l3-l3) + gl3*(6.*l1*l3-l1))*9./2.;
          dN[14] = (gl3*(3.*l2*l2-l2) + gl2*(6.*l2*l3-l3))*9./2.;
          dN[15] = (gl2*(3.*l3*l3-l3) + gl3*(6.*l2*l3-l2))*9./2.;
          dN[16] = (gl0*l1*l2 + gl1*l0*l2 + gl2*l0*l1)*27.;
          dN[17] = (gl0*l1*l3 + gl1*l0*l3 + gl3*l0*l1)*27.;
          dN[18] = (gl1*l2*l3 + gl2*l1*l3 + gl3*l1*l2)*27.;
          dN[19] = (gl0*l2*l3 + gl2*l0*l3 + gl3*l0*l2)*27.;
        }
        int countNodes() const {return 20;}
        void alignSharedNodes(Mesh* m, MeshEntity* elem,
            MeshEntity* edge, int order[]) {
          int which, rotate;
          bool flip;
          getAlignment(m, elem, edge, which, flip, rotate);
          if (! flip) {
            order[0] = 0;
            order[1] = 1;
          }
          else {
            order[0] = 1;
            order[1] = 0;
          }
        }
    };

    EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
      static Edge edge;
      static Triangle tri;
      static Tetrahedron tet;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,   // vertex
       &edge,     // edge
       &tri,      // triangle
       NULL,      // quad
       &tet,      // tet
       NULL,      // hex
       NULL,      // prism
       NULL};     // pyramid
      return shapes[type];
    }
    bool hasNodesIn(int dimension)
    {
      if (dimension < 3)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
        return 1;
      else if (type == Mesh::EDGE)
        return 2;
      else if (type == Mesh::TRIANGLE)
        return 1;
      else
        return 0;
    }
    int getOrder() {return 3;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
      assert(node < 2);
      if (type == Mesh::EDGE && node == 0)
        xi = Vector3(-1./3., 0, 0);
      else if (type == Mesh::EDGE && node == 1)
        xi = Vector3(1./3., 0, 0);
      else if (type == Mesh::TRIANGLE)
        xi = Vector3(1./3., 1./3., 0);
      else
        xi = Vector3(0, 0, 0);
    }
};

FieldShape* getLagrange(int order)
{
  static Linear linear;
  static LagrangeQuadratic quadratic;
  static LagrangeCubic cubic;
  if (order == 1)
    return &linear;
  if (order == 2)
    return &quadratic;
  if (order == 3)
    return &cubic;
  return NULL;
}

FieldShape* getSerendipity()
{
  static SerendipityQuadratic s;
  return &s;
}

/* these are step-wise fields which are defined by nodes
   at element centers, and the field value is constant
   throughout an element and discontinuous between elements.
   The first example is the gradient computed from a 1st-order
   Lagrange field, another example is the per-element error estimate */
template <int D>
class Constant : public FieldShape
{
  public:
    Constant()
    {
      std::stringstream ss;
      ss << "Constant_" << D;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const
    {
      return name.c_str();
    }
    class Element : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          grads.allocate(1);
          grads[0] = Vector3( 0, 0, 0);
        }
        int countNodes() const {return 1;}
        int getDimension() const {return D;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Element element;
      if (countNodesOn(type))
        return &element;
      return NULL;
    }
    bool hasNodesIn(int dimension)
    {
      if (dimension == D)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      int dimension = Mesh::typeDimension[type];
      if (dimension == D)
        return 1;
      else
        return 0;
   }
    int getOrder() {return 0;}
  private:
    std::string name;
};

FieldShape* getConstant(int dimension)
{
  static Constant<0> c0;
  static Constant<1> c1;
  static Constant<2> c2;
  static Constant<3> c3;
  static FieldShape* const table[4] =
  {&c0, &c1 ,&c2, &c3};
  return table[dimension];
}

int countElementNodes(FieldShape* s, int type)
{
  return s->getEntityShape(type)->countNodes();
}

}//namespace apf
