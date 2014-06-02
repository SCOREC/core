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

class Linear : public FieldShape
{
  public:
    const char* getName() const { return "Linear"; }
    class Vertex : public EntityShape
    {
      public:
        void getValues(Vector3 const& xi, NewArray<double>& values) const
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
      public:
        void getValues(Vector3 const&, NewArray<double>& ) const
        {
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>& ) const
        {
        }
        int countNodes() const {return 6;}
    };
    class Pyramid : public EntityShape
    {
      public:
        void getValues(Vector3 const&, NewArray<double>& ) const
        {
        }
        void getLocalGradients(Vector3 const&, NewArray<Vector3>& ) const
        {
        }
        int countNodes() const {return 5;}
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
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,      //vertex
       &edge,
       &triangle,
       &quad,
       &tet,
       NULL,      //hex
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

class Quadratic : public FieldShape
{
  public:
    const char* getName() const {return "Quadratic";}
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
      {&vertex,      //vertex
       &edge,     //edge
       &triangle, //triangle
       NULL,      //quad
       &tet,      //tet
       NULL,      //hex
       NULL,      //prism
       NULL};     //pyramid
      return shapes[type];
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
    int getOrder() {return 2;}
    void getNodeXi(int, int, Vector3& xi)
    {
      /* both for vertex nodes and mid-edge nodes, the xi
         coordinate is zero */
      xi = Vector3(0,0,0);
    }
};

FieldShape* getLagrange(int order)
{
  static Linear linear;
  static Quadratic quadratic;
  if (order == 1)
    return &linear;
  if (order == 2)
    return &quadratic;
  return NULL;
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

class IPShape : public FieldShape
{
  public:
    const char* getName() const {return name.c_str();}
    IPShape(int d, int o):dimension(d),order(o) {
      std::stringstream ss;
      ss << "IPShape_" << d << "_" << o;
      name = ss.str();
    }
    EntityShape* getEntityShape(int) {return 0;}
    bool hasNodesIn(int d)
    {
      return dimension==d;
    }
    int countNodesOn(int type)
    {
      if (Mesh::typeDimension[type]==dimension)
      {
        EntityIntegration const* ei = getIntegration(type);
        if (!ei) return 0; //currently some types don't exist, such as hex/prism
        return ei->getAccurate(order)->countPoints();
      }
      return 0;
    }
    /* the polynomial order of a set of integration points
       has no meaning */
    int getOrder() {return -1;}
  protected:
    int dimension;
    int order;
    std::string name;
};

FieldShape* getIPShape(int dimension, int order)
{
  static IPShape d3o1(3,1);
  static IPShape d3o2(3,2);
  static IPShape d3o3(3,3);
  static IPShape* table[4][4] =
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

class VoronoiShape : public IPShape
{
  public:
    const char* getName() const {return name.c_str();}
    VoronoiShape(int d, int o) :
      IPShape(d,o)
    {
      std::stringstream ss;
      ss << "IPShape_" << d << "_" << o;
      name = ss.str();
      elem.init(d,o);
    }
    EntityShape* getEntityShape(int type)
    {
      if (type == Mesh::TET)
        return &elem;
      return 0;
    }
    class Element : public EntityShape
    {
      public:
        void init(int d, int o)
        {
          dimension = d;
          order = o;
          /* BNG: Mesh:TET here and in CountNodes() should probably
             be generalized */
          EntityIntegration const* ei = getIntegration(Mesh::TET);
          int np = ei->getAccurate(order)->countPoints();
          points.setSize(np);
          for (int i = 0; i < np; ++i)
            points[i] = ei->getAccurate(order)->getPoint(i)->param;
        }
        int getClosestPtIdx(
            Vector3 const& p,
            DynamicArray<Vector3> const& points) const
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
        void getValues(
            Vector3 const& xi, 
            NewArray<double>& values) const
        {
          values.allocate(points.getSize());
          for (size_t i = 0; i < points.getSize(); ++i)
            values[i] = 0.0;
          values[getClosestPtIdx(xi,points)] = 1.0;
        }
        void getLocalGradients(
            Vector3 const& xi, 
            NewArray<Vector3>& grads) const
        {
          fail("gradients not defined for Voronoi shapes");
        }
        int countNodes() const
        {
          return points.getSize();
        }
        int dimension;
        int order;
        std::string name;
        DynamicArray<Vector3> points;
    };
    void getNodeXi(int type, int node, Vector3& xi)
    {
      assert(type == Mesh::TET);
      xi = elem.points[node];
    }
  private:
    Element elem;
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
