/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include <vector>

namespace apf {

static const double c = -2.44948974278318;

class Hierarchic : public FieldShape
{
  public:
    Hierarchic() {}
    const char* getName() const { return "Hierarchic"; }
    class Vertex : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& N) const
        {
          N.allocate(1);
          N[0] = 1.0;
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
            Vector3 const& xi, NewArray<double>& N) const
        {
          N.allocate(3);
          N[0] = (1.0-xi[0])/2.0;
          N[1] = (1.0+xi[0])/2.0;
          N[2] = c*N[0]*N[1];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          dN.allocate(3);
          dN[0] = Vector3(-0.5, 0.0, 0.0);
          dN[1] = Vector3( 0.5, 0.0, 0.0);
          dN[2] = Vector3(-0.5*c*xi[0], 0.0, 0.0);
        }
        int countNodes() const {return 3;}
    };
    class Triangle : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& N) const
        {
          N.allocate(6);
          N[0] = 1.0-xi[0]-xi[1];
          N[1] = xi[0];
          N[2] = xi[1];
          N[3] = c*N[0]*N[1];
          N[4] = c*N[1]*N[2];
          N[5] = c*N[2]*N[0];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          dN.allocate(6);
          dN[0] = Vector3(-1.0, -1.0, 0.0);
          dN[1] = Vector3( 1.0,  0.0, 0.0);
          dN[2] = Vector3( 0.0,  1.0, 0.0);
          dN[3] = Vector3( 1.0-2.0*xi[0]-xi[1],  -xi[0], 0.0) * c;
          dN[4] = Vector3( xi[1], xi[0], 0.0) * c;
          dN[5] = Vector3( -xi[1],  1.0-xi[0]-2.0*xi[1], 0.0) * c;
        }
        int countNodes() const {return 6;}
    };
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& N) const
        {
          N.allocate(10);
          N[0] = 1.0-xi[0]-xi[1]-xi[2];
          N[1] = xi[0];
          N[2] = xi[1];
          N[3] = xi[2];
          N[4] = c*N[0]*N[1];
          N[5] = c*N[1]*N[2];
          N[6] = c*N[2]*N[0];
          N[7] = c*N[0]*N[3];
          N[8] = c*N[1]*N[3];
          N[9] = c*N[2]*N[3];
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          dN.allocate(10);
          dN[0] = Vector3(-1.0, -1.0, -1.0);
          dN[1] = Vector3( 1.0,  0.0,  0.0);
          dN[2] = Vector3( 0.0,  1.0,  0.0);
          dN[3] = Vector3( 0.0,  0.0,  1.0);
          dN[4] = Vector3( 1.0-2.0*xi[0]-xi[1]-xi[2], -xi[0], -xi[0] ) * c;
          dN[5] = Vector3( xi[1], xi[0], 0.0 ) * c;
          dN[6] = Vector3( -xi[1], 1.0-xi[0]-2.0*xi[1]-xi[2], -xi[1] ) * c;
          dN[7] = Vector3( -xi[2], -xi[2], 1.0-xi[0]-xi[1]-2.0*xi[2] ) * c;
          dN[8] = Vector3( xi[2], 0.0, xi[0] ) * c;
          dN[9] = Vector3( 0.0, xi[2], xi[1] ) * c;
        }
        int countNodes() const {return 10;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
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
    void getNodeXi(int, int, Vector3&)
    {
      /* this has no meaning for higher order hierarchic nodes.
         My guess is ma uses it for solution transfer, for which
         something new would need to happen for this shape. */
      fail("unimplemented getNodeXi called");
    }
    bool hasNodesIn(int dimension)
    {
      if ( (dimension == 0) || (dimension == 1) )
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if ( (type == Mesh::VERTEX) || (type == Mesh::EDGE) )
        return 1;
      else
        return 0;
    }
    int getOrder() {return 2;}
};

FieldShape* getHierarchic(int o)
{
  static Hierarchic q;
  if (o == 1)
    return getLagrange(o);
  else if (o == 2)
    return &q;
  return NULL;
}

}
