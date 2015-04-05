/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfIntegrate.h"

namespace apf {

static std::map<std::string, FieldShape*> registry;

EntityShape::~EntityShape()
{
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
        void getValues(Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1.0;
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>&) const
        {
        }
        int countNodes() const {return 1;}
    };
    class Edge : public EntityShape
    {
      public:
        void getValues(Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(2);
          values[0] = (1.0-xi[0])/2.0;
          values[1] = (1.0+xi[0])/2.0;
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>& grads) const
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(3);
          values[0] = 1-xi[0]-xi[1];
          values[1] = xi[0];
          values[2] = xi[1];
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>& grads) const
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
        void getLocalGradients(Vector3 const& xi,
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(4);
          values[0] = 1-xi[0]-xi[1]-xi[2];
          values[1] = xi[0];
          values[2] = xi[1];
          values[3] = xi[2];
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>& grads) const
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
        void getLocalGradients(Vector3 const& xi,
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
        void getLocalGradients(Vector3 const& xi,
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
          values[4] = l0x * l0y * l1z / 8;
          values[5] = l1x * l0y * l1z / 8;
          values[6] = l1x * l1y * l1z / 8;
          values[7] = l0x * l1y * l1z / 8;
        }
        void getLocalGradients(Vector3 const& xi,
            NewArray<Vector3>& grads) const
        {
          double l0x = (1 - xi[0]);
          double l1x = (1 + xi[0]);
          double l0y = (1 - xi[1]);
          double l1y = (1 + xi[1]);
          double l0z = (1 - xi[2]);
          double l1z = (1 + xi[2]);
          grads.allocate(5);
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
};

class QuadraticBase : public FieldShape
{
  public:
    class Edge : public EntityShape
    {
      public:
        void getValues(Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(3);
          values[0] = -xi[0]*(1-xi[0])/2.0;
          values[1] =  xi[0]*(1+xi[0])/2.0;
          values[2] = 1-xi[0]*xi[0];
        }
        void getLocalGradients(Vector3 const& xi,
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
        void getLocalGradients(Vector3 const& xi,
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
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
        void getLocalGradients(Vector3 const& xi,
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
    { /* TODO: implement this and then update hasNodesIn, countNodesOn */
      public:
        void getValues(Vector3 const&, NewArray<double>&) const
        {
          fail("quadratic Lagrange quadrilateral shape values not implemented\n");
        }
        void getLocalGradients(Vector3 const&,
            NewArray<Vector3>&) const
        {
          fail("quadratic Lagrange quadrilateral shape grads not implemented\n");
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

class SerendipityQuadratic : public QuadraticBase
{
  public:
    SerendipityQuadratic() { registerSelf(apf::SerendipityQuadratic::getName()); }
    const char* getName() const {return "Serendipity Quadratic";}
    class Quad : public EntityShape
    {
    public:
      void getValues(Vector3 const& xi, NewArray<double>& values) const
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
      void getLocalGradients(Vector3 const& xi,
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

FieldShape* getLagrange(int order)
{
  static Linear linear;
  static LagrangeQuadratic quadratic;
  if (order == 1)
    return &linear;
  if (order == 2)
    return &quadratic;
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
        void getValues(Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1;
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>& grads) const
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

static Integration const* tryToGetIntegration(int type, int order)
{
  EntityIntegration const* ei = getIntegration(type);
  if (!ei)
    return 0;
  return ei->getAccurate(order);
}

class IPBase : public FieldShape
{
  public:
    IPBase(int d, int o):dimension(d),order(o) {
    }
    EntityShape* getEntityShape(int) {return 0;}
    bool hasNodesIn(int d)
    {
      return dimension==d;
    }
    int countNodesOn(int type)
    {
      if (Mesh::typeDimension[type]!=dimension)
        return 0; //non-elements have no integration points
      Integration const* i = tryToGetIntegration(type,order);
      if (!i)
        return 0; //some types have no integrations of this order
      return i->countPoints();
    }
    /* this field can integrate a polynomial of
       this order exactly */
    int getOrder() {return order;}
  protected:
    int dimension;
    int order;
};

class IPShape : public IPBase
{
  public:
    IPShape(int d, int o):
      IPBase(d, o)
    {
      std::stringstream ss;
      ss << "IPShape_" << d << "_" << o;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const {return name.c_str();}
  private:
    std::string name;
};

FieldShape* getIPShape(int dimension, int order)
{
  static IPShape d2o1(2,1);
  static IPShape d2o2(2,2);
  static IPShape d2o3(2,3);
  static IPShape d3o1(3,1);
  static IPShape d3o2(3,2);
  static IPShape d3o3(3,3);
  static IPShape* table[4][4] =
  {{0,0,0,0}//vertex
  ,{0,0,0,0}//edge
  ,{0,&d2o1,&d2o2,&d2o3}//face
  ,{0,&d3o1,&d3o2,&d3o3}//region
  };
  assert(dimension >= 0);
  assert(dimension <= 3);
  assert(order >= 0);
  assert(order <= 3);
  return table[dimension][order];
}

class VoronoiShape : public IPBase
{
  public:
    VoronoiShape(int d, int o) :
      IPBase(d,o)
    {
      std::stringstream ss;
      ss << "VoronoiShape_" << d << "_" << o;
      name = ss.str();
      registerSelf(name.c_str());
      for (int type = 0; type < Mesh::TYPES; ++type)
        if (Mesh::typeDimension[type] == d)
          elem[type].init(type,o);
    }
    const char* getName() const {return name.c_str();}
    EntityShape* getEntityShape(int type)
    {
      return &elem[type];
    }
    class Element : public EntityShape
    {
      public:
        void init(int type, int order)
        {
          Integration const* in = tryToGetIntegration(type,order);
          if (!in)
            return;
          int np = in->countPoints();
          points.setSize(np);
          for (int i = 0; i < np; ++i)
            points[i] = in->getPoint(i)->param;
        }
        int getClosestPtIdx(Vector3 const& p) const
        {
          int idx = 0;
          double leastDistance = (p - points[0]).getLength();
          for (size_t i = 1; i < points.getSize(); ++i)
          {
            double distance = (p - points[i]).getLength();
            if (distance < leastDistance)
            {
              leastDistance = distance;
              idx = i;
            }
          }
          return idx;
        }
        void getValues(Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(points.getSize());
          for (size_t i = 0; i < points.getSize(); ++i)
            values[i] = 0.0;
          values[getClosestPtIdx(xi)] = 1.0;
        }
        void getLocalGradients(
            Vector3 const&,
            NewArray<Vector3>&) const
        {
          fail("gradients not defined for Voronoi shapes");
        }
        int countNodes() const
        {
          return points.getSize();
        }
        DynamicArray<Vector3> points;
    };
    void getNodeXi(int type, int node, Vector3& xi)
    {
      xi = elem[type].points[node];
    }
  private:
    Element elem[Mesh::TYPES];
    std::string name;
};
 
FieldShape* getVoronoiShape(int dimension, int order)
{
  static VoronoiShape d3o1(3,1);
  static VoronoiShape d3o2(3,2);
  static VoronoiShape d3o3(3,3);
  static VoronoiShape* table[4][4] =
  {{0,0,0,0}//vertex
  ,{0,0,0,0}//edge
  ,{0,0,0,0}//face
  ,{0,&d3o1,&d3o2,&d3o3}//region
  };
  assert(dimension >= 0);
  assert(dimension <= 3);
  assert(order >= 0);
  assert(order <= 3);
  return table[dimension][order];
}

int countElementNodes(FieldShape* s, int type)
{
  return s->getEntityShape(type)->countNodes();
}

}//namespace apf
