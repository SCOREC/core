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

class Mode
{
  public:
    virtual ~Mode() {};
    virtual void getVals(Vector3 const& xi, NewArray<double>& N) const = 0;
    virtual void getGrads(Vector3 const& xi, NewArray<Vector3>& dN) const = 0;
    virtual int countNodes() const = 0;
    void pushVals(Vector3 const& xi, std::vector<double>& N)
    {
      NewArray<double> NN;
      this->getVals(xi,NN);
      for (int i=0; i < this->countNodes(); ++i)
        N.push_back(NN[i]);
    }
    void pushGrads(Vector3 const& xi, std::vector<Vector3>& dN)
    {
      NewArray<Vector3> dNN;
      this->getGrads(xi,dNN);
      for (int i=0; i < this->countNodes(); ++i)
        dN.push_back(dNN[i]);
    }
};

class VertexMode : public Mode
{
  public:
    void getVals(Vector3 const& xi, NewArray<double>& N) const
    {
      return getTetCoord(xi,N);
    }
    void getGrads(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
      return getTetCoordGrad(dN);
    }
    int countNodes() const
    {
      return 4;
    }
};

class EdgeMode : public Mode
{
  public:
    EdgeMode(int o) : p(o)
    {
      assert( p>=2 );
    }
    void getVals(Vector3 const& xi, NewArray<double>& N) const
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
    void getGrads(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
      NewArray<double> l;
      getTetCoord(xi,l);
      NewArray<Vector3> dl;
      getTetCoordGrad(dl);
      dN.allocate(6);
      for (int i=0; i < 3; ++i)
      {
        dN[0][i] = phi(p-2,l[0] -l[1])  * (l[0] * dl[1][i] + dl[0][i] * l[1]) +
          dphi(p-2,l[0] - l[1]) * (dl[0][i] - dl[1][i]);
        dN[1][i] = phi(p-2,l[1] -l[2])  * (l[1] * dl[2][i] + dl[1][i] * l[2]) +
          dphi(p-2,l[1] - l[2]) * (dl[1][i] - dl[2][i]);
        dN[2][i] = phi(p-2,l[2] -l[0])  * (l[2] * dl[0][i] + dl[2][i] * l[0]) +
          dphi(p-2,l[2] - l[0]) * (dl[2][i] - dl[0][i]);
        dN[3][i] = phi(p-2,l[0] -l[3])  * (l[0] * dl[3][i] + dl[0][i] * l[3]) +
          dphi(p-2,l[0] - l[3]) * (dl[0][i] - dl[3][i]);
        dN[4][i] = phi(p-2,l[1] -l[3])  * (l[1] * dl[3][i] + dl[1][i] * l[3]) +
          dphi(p-2,l[1] - l[3]) * (dl[1][i] - dl[3][i]);
        dN[5][i] = phi(p-2,l[2] -l[3])  * (l[2] * dl[3][i] + dl[2][i] * l[3]) +
          dphi(p-2,l[2] - l[3]) * (dl[2][i] - dl[3][i]);
      }
    }
    int countNodes() const
    {
      return 6;
    }
  private:
    int p;
};

class FaceMode : public Mode
{
  public:
    FaceMode(int o) : p(o)
    {
      assert (p >= 3);
      fprintf(stderr,"unimplemented face mode");
    }
    void getVals(Vector3 const& xi, NewArray<double>& N) const
    {
    }
    void getGrads(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
    }
    int countNodes() const
    {
      return 4 * (p-2);
    }
  private:
    int p;
};

class RegionMode : public Mode
{
  public:
    RegionMode(int o) : p(o)
    {
      assert (p >= 4);
      fprintf(stderr,"unimplemented region mode");
    }
    void getVals(Vector3 const& xi, NewArray<double>& N) const
    {
    }
    void getGrads(Vector3 const& xi, NewArray<Vector3>& dN) const
    {
    }
    int countNodes() const
    {
      return (p-2) * (p-3) / 2;
    }
  private:
    int p;
};

class Hierarchic : public FieldShape
{
  public:
    Hierarchic(int order) : p(order)
    {
      assert (p >=1 );
      assert (p <= 2);
    }
    const char* getName() const { return "Hierarchic"; }
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
    class Tetrahedron : public EntityShape
    {
      public:
        Tetrahedron(int o) : p(o) {}
        void getValues(Vector3 const& xi, NewArray<double>& N) const
        {
          std::vector<double> NN;
          VertexMode v;
          v.pushVals(xi,NN);
          for (int o=2; o <= p; ++o)
          {
            EdgeMode e(o);
            e.pushVals(xi,NN);
            if (o >= 3)
            {
              FaceMode f(o);
              f.pushVals(xi,NN);
            }
            if (o >= 4)
            {
              RegionMode r(o);
              r.pushVals(xi,NN);
            }
          }
          int n = this->countNodes();
          assert(NN.size() == n);
          N.allocate(n);
          for (int i=0; i < n; ++i)
            N[i] = NN[i];
        }
        void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& dN) const
        {
          std::vector<Vector3> dNN;
          VertexMode v;
          v.pushGrads(xi,dNN);
          for (int o=2; o <= p; ++o)
          {
            EdgeMode e(o);
            e.pushGrads(xi,dNN);
            if (o >= 3)
            {
              FaceMode f(o);
              f.pushGrads(xi,dNN);
            }
            if (o >= 4)
            {
              RegionMode r(o);
              r.pushGrads(xi,dNN);
            }
          }
          int n = this->countNodes();
          assert(dNN.size() == n);
          dN.allocate(n);
          for (int i=0; i < n; ++i)
            dN[i] = dNN[i];
        }
        int countNodes() const
        {
          VertexMode v;
          int n = v.countNodes();
          for (int o=2; o <= p; ++o)
          {
            EdgeMode e(o);
            n += e.countNodes();
            if (o >=3)
            {
              FaceMode f(o);
              n += f.countNodes();
            }
            if (o >= 4)
            {
              RegionMode r(o);
              n += r.countNodes();
            }
          }
          return n;
        }
      private:
        int p;
    };
    EntityShape* getEntityShape(int type)
    {
      static Vertex vtx;
      static Tetrahedron tet(p);
      static EntityShape* shapes[Mesh::TYPES] =
      {&vtx,  //vertex
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
      else if ( (dimension == 1) && (p >= 2) )
        return true;
      else if ( (dimension == 2) && (p >= 3) )
        return true;
      else if ( (dimension == 3) && (p >= 4) )
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
        return 1;
      else if (type == Mesh::EDGE)
        return p-1;
      else if (type == Mesh::TRIANGLE)
        return (p-1) * (p-2) / 2;
      else if (type == Mesh::TET)
        return (p-1) * (p-2) * (p-3) / 6;
      else
        return 0;
    }
    int getOrder() { return p; }
    void getNodeXi(int, int, Vector3& xi)
    {
      xi = Vector3(0,0,0);
    }
  private:
    int p;
};

FieldShape* getHierarchic(int order)
{
  assert( (order == 1) || (order == 2) );
  static Hierarchic hierarchic(order);
  return &hierarchic;
}

apf::Field* createHierarchicField(apf::Mesh* m, const char* name,
    int valueType, int order)
{
  return createField(m,name,valueType,getHierarchic(order));
}

}
