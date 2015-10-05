/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfIntegrate.h"
#include <canArray.h>
#include <cassert>

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
  ,{0,&d2o1,&d2o2,&d2o3,&d2o4,0,0,0}//face
  ,{0,&d3o1,&d3o2,&d3o3,&d3o4,&d3o5,&d3o6,&d3o7}//region
  };
  assert(dimension >= 0);
  assert(dimension <= 3);
  assert(order >= 0);
  assert(order <= 7);
  return table[dimension][order];
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
  static VoronoiShape d2o1(2,1);
  static VoronoiShape d2o2(2,2);
  static VoronoiShape d2o3(2,3);
  static VoronoiShape d3o1(3,1);
  static VoronoiShape d3o2(3,2);
  static VoronoiShape d3o3(3,3);
  static VoronoiShape* table[4][4] =
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

}
