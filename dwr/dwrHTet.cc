/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>

namespace dwr {

using namespace apf;

static const double c0 = -2.44948974278318;  /* -2 sqrt(3/2) */
static const double c1 = -3.16227766016838;  /* -2 sqrt(5/2) */
static const double c2 = -0.935414346693485; /* -1/2 sqrt(7/2) */
static const double c3 = -1.06066017177982;  /* -1/2 sqrt(9/2) */

static double phi(int o, double x)
{
  double v = 0.0;
  switch (o)
  {
    case 0:
      v = c0;
      break;
    case 1:
      v = c1 * x;
      break;
    case 2:
      v = c2 * (5.0*x*x - 1.0);
      break;
    case 3:
      v = c3 * (7.0*x*x - 3.0) * x;
      break;
    default:
      fprintf(stderr,"unsupported order");
  }
  return v;
}

static double dphi(int o, double x)
{
  double v = 0.0;
  switch (o)
  {
    case 0:
      v = 0.0;
      break;
    case 1:
      v = c1;
      break;
    case 2:
      v = c2 * 10.0 * x;
      break;
    case 3:
      v = c3 * (21.0*x*x - 3.0);
      break;
    default:
      fprintf(stderr,"unsupoorted order");
  }
  return v;
}

static void getTetCoord(Vector3 const& xi, NewArray<double>& l)
{
  l.allocate(4);
  l[0] = -0.5 * (xi[0] + xi[1] + xi[2] + 1.0);
  l[1] = 0.5 * (xi[0] + 1.0);
  l[2] = 0.5 * (xi[1] + 1.0);
  l[3] = 0.5 * (xi[2] + 1.0);
}

static void getTetCoordGrad(NewArray<Vector3>& dl)
{
  dl.allocate(4);
  dl[0] = Vector3(-0.5, -0.5, -0.5);
  dl[1] = Vector3( 0.5,  0.0,  0.0);
  dl[2] = Vector3( 0.0,  0.5,  0.0);
  dl[3] = Vector3( 0.0,  0.0,  0.5);
}

class VertexMode : public EntityShape
{
  public:
    void getValues(Vector3 const& xi, NewArray<double>& N) const
    {
      getTetCoord(xi,N);
    }
    void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
      getTetCoordGrad(dN);
    }
    int countNodes() const
    {
      return 4;
    }
};

class EdgeMode : public EntityShape
{
  public:
    EdgeMode(int o) : p(o)
    {
      assert( p>=2 );
    }
    void getValues(Vector3 const& xi, NewArray<double>& N) const
    {
      NewArray<double> l;
      getTetCoord(xi,l);
      N.allocate(6);
      N[0] = l[0] * l[1] * phi(p-2, l[0] - l[1]);
      N[1] = l[1] * l[2] * phi(p-2, l[1] - l[2]);
      N[2] = l[2] * l[0] * phi(p-2, l[2] - l[0]);
      N[3] = l[0] * l[3] * phi(p-2, l[0] - l[3]);
      N[4] = l[1] * l[3] * phi(p-2, l[1] - l[3]);
      N[5] = l[2] * l[3] * phi(p-2, l[2] - l[3]);
    }
    void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
      NewArray<double> l;
      getTetCoord(xi,l);
      NewArray<Vector3> dl;
      getTetCoordGrad(dl);
      dN.allocate(6);
      for (int i=0; i < 3; ++i)
      {
        dN[0][i] = phi(p-2, l[0] -l[1])  * (l[0] * dl[1][i] + dl[0][i] * l[1]) +
          dphi(p-2, l[0] - l[1]) * (dl[0][i] - dl[1][i]);
        dN[1][i] = phi(p-2, l[1] -l[2])  * (l[1] * dl[2][i] + dl[1][i] * l[2]) +
          dphi(p-2, l[1] - l[2]) * (dl[1][i] - dl[2][i]);
        dN[2][i] = phi(p-2, l[2] -l[0])  * (l[2] * dl[0][i] + dl[2][i] * l[0]) +
          dphi(p-2, l[2] - l[0]) * (dl[2][i] - dl[0][i]);
        dN[3][i] = phi(p-2, l[0] -l[3])  * (l[0] * dl[3][i] + dl[0][i] * l[3]) +
          dphi(p-2, l[0] - l[3]) * (dl[0][i] - dl[3][i]);
        dN[4][i] = phi(p-2, l[1] -l[3])  * (l[1] * dl[3][i] + dl[1][i] * l[3]) +
          dphi(p-2, l[1] - l[3]) * (dl[1][i] - dl[3][i]);
        dN[5][i] = phi(p-2, l[2] -l[3])  * (l[2] * dl[3][i] + dl[2][i] * l[3]) +
          dphi(p-2, l[2] - l[3]) * (dl[2][i] - dl[3][i]);
      }
    }
    int countNodes() const
    {
      return 6;
    }
  private:
    int p;
};

class FaceMode : public EntityShape
{
  public:
    FaceMode(int o) : p(o)
    {
      assert (p >= 3);
      fprintf(stderr,"unimplemented face mode");
    }
    void getValues(Vector3 const& xi, NewArray<double>& N) const
    {
    }
    void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
    }
    int countNodes() const
    {
      return 4 * (p-2);
    }
  private:
    int p;
};

class RegionMode : public EntityShape
{
  public:
    RegionMode(int o) : p(o)
    {
      assert (p >= 4);
      fprintf(stderr,"unimplemented region mode");
    }
    void getValues(Vector3 const& xi, NewArray<double>& N) const
    {
    }
    void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
    }
    int countNodes() const
    {
      return (p-2) * (p-3) / 2;
    }
  private:
    int p;
};

class HLinear : public FieldShape
{
  public:
    const char* getName() const { return "HLinear"; }
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Vector3 const& xi, NewArray<double>& N) const
        {
          VertexMode v;
          v.getValues(xi,N);
        }
        void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          VertexMode v;
          v.getLocalGradients(xi,dN);
        }
        int countNodes() const
        {
          VertexMode v;
          return v.countNodes();
        }
    };
    EntityShape* getEntityShape(int type)
    {
      static Tetrahedron tet;
      static EntityShape* shapes[Mesh::TYPES] =
      {NULL,  //vertex
       NULL,  //edge
       NULL,  //triangle
       NULL,  //quad
       &tet,  //tet
       NULL,  //hex
       NULL,  //prism
       NULL}; //pyramid
      return shapes[type];
    }
    bool hasNodesIn(int dimension)
    {
      if ( dimension == 0 )
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
    int getOrder() { return 2; }
    /* find out what this means and why it exists*/
    void getNodeXi(int, int, Vector3& xi)
    {
      xi = Vector3(0,0,0);
    }
};

FieldShape* getHierarchic(int order)
{
  static HLinear linear;
  if (order == 1)
    return &linear;
  return NULL;
}

}
