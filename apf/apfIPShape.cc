/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfIntegrate.h"
#include <mthMatrix.h>
#include <pcu_util.h>

#include <iostream>

namespace apf {

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
  static IPShape d2o4(2,4);
  static IPShape d2o5(2,5);
  static IPShape d3o1(3,1);
  static IPShape d3o2(3,2);
  static IPShape d3o3(3,3);
  static IPShape d3o4(3,4);
  static IPShape d3o5(3,5);
  static IPShape d3o6(3,6);
  static IPShape d3o7(3,7);
  static IPShape* table[4][8] =
  {{0,0,0,0,0}//vertex
  ,{0,0,0,0,0}//edge
  ,{0,&d2o1,&d2o2,&d2o3,&d2o4,&d2o5,0,0}//face
  ,{0,&d3o1,&d3o2,&d3o3,&d3o4,&d3o5,&d3o6,&d3o7}//region
  };
  PCU_ALWAYS_ASSERT(dimension >= 0);
  PCU_ALWAYS_ASSERT(dimension <= 3);
  PCU_ALWAYS_ASSERT(order >= 0);
  PCU_ALWAYS_ASSERT(order <= 7);
  apf::FieldShape* shape = table[dimension][order];
  PCU_ALWAYS_ASSERT(shape);
  return shape;
}

static void getIntegrationPoints(
    int type, int order, can::Array<Vector3>& points)
{
  Integration const* in = tryToGetIntegration(type,order);
  if (!in)
    return;
  int np = in->countPoints();
  points.resize(np);
  for (int i = 0; i < np; ++i)
    points[i] = in->getPoint(i)->param;
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
          getIntegrationPoints(type, order, points);
        }
        int getClosestPtIdx(Vector3 const& p) const
        {
          int idx = 0;
          double leastDistance = (p - points[0]).getLength();
          for (size_t i = 1; i < points.size(); ++i)
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
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(points.size());
          for (size_t i = 0; i < points.size(); ++i)
            values[i] = 0.0;
          values[getClosestPtIdx(xi)] = 1.0;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&,
            NewArray<Vector3>&) const
        {
          fail("gradients not defined for Voronoi shapes");
        }
        int countNodes() const
        {
          return points.size();
        }
        can::Array<Vector3> points;
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
  static VoronoiShape d1o1(1,1);
  static VoronoiShape d2o1(2,1);
  static VoronoiShape d2o2(2,2);
  static VoronoiShape d2o3(2,3);
  static VoronoiShape d2o4(2,4);
  static VoronoiShape d3o1(3,1);
  static VoronoiShape d3o2(3,2);
  static VoronoiShape d3o3(3,3);
  static VoronoiShape d3o4(3,4);
  static VoronoiShape* table[4][5] =
  {{0,0,0,0,0}//vertex
  ,{0,&d1o1,0,0,0}//edge
  ,{0,&d2o1,&d2o2,&d2o3,&d2o4}//face
  ,{0,&d3o1,&d3o2,&d3o3,&d3o4}//region
  };
  PCU_ALWAYS_ASSERT(dimension >= 0);
  PCU_ALWAYS_ASSERT(dimension <= 3);
  PCU_ALWAYS_ASSERT(order >= 0);
  PCU_ALWAYS_ASSERT(order <= 4);
  return table[dimension][order];
}

class ConstantIPFit : public IPBase
{
  public:
    ConstantIPFit(int d) : IPBase(d, 1)
    {
      std::stringstream ss;
      ss << "ConstantIPFit_" << d;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const {return name.c_str();}
    class Triangle : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1.0;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          grads.allocate(1);
          grads[0] = Vector3(0,0,0);
        }
        int countNodes() const {return 1;}
    };
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
          values.allocate(1);
          values[0] = 1.0;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          grads.allocate(1);
          grads[0] = Vector3(0,0,0);
        }
        int countNodes() const {return 1;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Triangle triangle;
      static Tetrahedron tet;
      static EntityShape* shapes[Mesh::TYPES] =
      {0,         // vertex
       0,         // edge
       &triangle, // triangle
       0,         // quad
       &tet,      // tet
       0,         // prism
       0,         // pyramid
       0};        // hex
      return shapes[type];
    }
    void getNodeXi(int type, int node, Vector3& xi)
    {
      can::Array<Vector3> points;
      getIntegrationPoints(type, 1, points);
      xi = points[node];
    }
  private:
    std::string name;
};

class LinearIPFit : public IPBase
{
  public:
    LinearIPFit(int d) : IPBase(d, 2)
    {
      std::stringstream ss;
      ss << "LinearIPFit_" << d;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const {return name.c_str();}
    class Triangle : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(3);

          mth::Matrix<double,3,3> c;
          c(0,0) = -0.333333333333334; c(0,1) = 2.000000000000000; c(0,2) = 0.000000000000000;
          c(1,0) = -0.333333333333334; c(1,1) = 0.000000000000000; c(1,2) = 2.000000000000000;
          c(2,0) = 1.666666666666668; c(2,1) = -2.000000000000000; c(2,2) = -2.000000000000000;

          mth::Vector<double,3> p;
          p(0) = 1.0;
          p(1) = xi[0];
          p(2) = xi[1];

          for (int i=0; i < 3; ++i)
            values[i] = c[i] * p;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
          fail("grads not implemented yet");
        }
        int countNodes() const {return 3;}
    };
    class Tetrahedron : public EntityShape
    {
      public:
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          values.allocate(4);

          mth::Matrix<double,4,4> c;
          c(0,0) = 1.927050983124845; c(0,1) = -2.236067977499789; c(0,2) = -2.236067977499789; c(0,3) = -2.236067977499789;
          c(1,0) = -0.309016994374948; c(1,1) = 2.236067977499789; c(1,2) = 0.000000000000000; c(1,3) = 0.000000000000000;
          c(2,0) = -0.309016994374948; c(2,1) = 0.000000000000000; c(2,2) = 2.236067977499789; c(2,3) = 0.000000000000000;
          c(3,0) = -0.309016994374948; c(3,1) = 0.000000000000000; c(3,2) = 0.000000000000000; c(3,3) = 2.236067977499789;

          mth::Vector<double,4> p;
          p(0) = 1.0;
          p(1) = xi[0];
          p(2) = xi[1];
          p(3) = xi[2];

          for (int i=0; i < 4; ++i)
            values[i] = c[i] * p;
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
          fail("grads not implemented yet");
        }
        int countNodes() const {return 4;}
    };
    EntityShape* getEntityShape(int type)
    {
      static Triangle triangle;
      static Tetrahedron tet;
      static EntityShape* shapes[Mesh::TYPES] =
      {0,         // vertex
       0,         // edge
       &triangle, // triangle
       0,         // quad
       &tet,      // tet
       0,         // prism
       0,         // pyramid
       0};        // hex
      return shapes[type];
    }
    void getNodeXi(int type, int node, Vector3& xi)
    {
      can::Array<Vector3> points;
      getIntegrationPoints(type, 2, points);
      xi = points[node];
    }
  private:
    std::string name;
};

FieldShape* getIPFitShape(int dimension, int order)
{
  static ConstantIPFit d2o1(2);
  static ConstantIPFit d3o1(3);
  static LinearIPFit d2o2(2);
  static LinearIPFit d3o2(3);
  FieldShape* table[4][3] =
  {{0,0,0}//vertex
  ,{0,0,0}//edge
  ,{0,&d2o1,&d2o2}//face
  ,{0,&d3o1,&d3o2}};//region
  PCU_ALWAYS_ASSERT(dimension >= 0);
  PCU_ALWAYS_ASSERT(dimension <= 3);
  PCU_ALWAYS_ASSERT(order >= 0);
  PCU_ALWAYS_ASSERT(order <= 2);
  return table[dimension][order];
}

}
