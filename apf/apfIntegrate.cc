/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfIntegrate.h"
#include "apfMesh.h"
#include "apf.h"

namespace apf {

Integration const* EntityIntegration::getAccurate(int minimumAccuracy) const
{
  int n = this->countIntegrations();
  for (int i=0; i < n; ++i)
  {
    Integration const* integration = this->getIntegration(i);
    if (integration->getAccuracy() >= minimumAccuracy)
      return integration;
  }
  return NULL;
}

class EdgeIntegration : public EntityIntegration
{
  public:
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          static IntegrationPoint point(Vector3(0,0,0),2);
          return &point;
        }
        virtual int getAccuracy() const {return 1;}
    };
    class N2 : public Integration
    {
      public:
        virtual int countPoints() const {return 2;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[2]=
          { IntegrationPoint(Vector3( 0.577350269189626,0,0),1),
            IntegrationPoint(Vector3(-0.577350269189626,0,0),1) };
          return points + i;
        }
        virtual int getAccuracy() const {return 3;}
    };
    class N3 : public Integration
    {
      public:
        virtual int countPoints() const {return 3;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[3]=
          { IntegrationPoint(Vector3( 0.000000000000000,0,0),0.888888888888889),
            IntegrationPoint(Vector3( 0.774596669241483,0,0),0.555555555555556),
            IntegrationPoint(Vector3(-0.774596669241483,0,0),0.555555555555556) };
          return points + i;
        }
        virtual int getAccuracy() const {return 5;}
    };
    class N4 : public Integration
    {
      public:
        virtual int countPoints() const {return 4;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[4]=
          { IntegrationPoint(Vector3( 0.339981043584856,0,0),0.652145154862546),
            IntegrationPoint(Vector3(-0.339981043584856,0,0),0.652145154862546),
            IntegrationPoint(Vector3( 0.861136311594053,0,0),0.347854845137454),
            IntegrationPoint(Vector3(-0.861136311594053,0,0),0.347854845137454) };
          return points + i;
        }
        virtual int getAccuracy() const {return 7;}
    };
    class N5 : public Integration
    {
      public:
        virtual int countPoints() const {return 5;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[5]=
          { IntegrationPoint(Vector3( 0.000000000000000,0,0),0.568888888888889),
            IntegrationPoint(Vector3( 0.538469310105683,0,0),0.478628670499366),
            IntegrationPoint(Vector3(-0.538469310105683,0,0),0.478628670499366),
            IntegrationPoint(Vector3( 0.906179845938664,0,0),0.236926885056189),
            IntegrationPoint(Vector3(-0.906179845938664,0,0),0.236926885056189) };
          return points + i;
        }
        virtual int getAccuracy() const {return 9;}
    };
    class N6 : public Integration
    {
      public:
        virtual int countPoints() const {return 6;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[6]=
          { IntegrationPoint(Vector3( 0.238619186083197,0,0),0.467913934572691),
            IntegrationPoint(Vector3(-0.238619186083197,0,0),0.467913934572691),
            IntegrationPoint(Vector3( 0.661209386466265,0,0),0.360761573048139),
            IntegrationPoint(Vector3(-0.661209386466265,0,0),0.360761573048139),
            IntegrationPoint(Vector3( 0.932469514203152,0,0),0.171324492379170),
            IntegrationPoint(Vector3(-0.932469514203152,0,0),0.171324492379170) };
          return points + i;
        }
        virtual int getAccuracy() const {return 11;}
    };
    virtual int countIntegrations() const {return 6;}
    virtual Integration const* getIntegration(int i) const
    {
      static N1 i1;
      static N2 i2;
      static N3 i3;
      static N4 i4;
      static N5 i5;
      static N6 i6;
      static Integration* integrations[6] = 
      {&i1,&i2,&i3,&i4,&i5,&i6};
      return integrations[i];
    }
};

class TriangleIntegration : public EntityIntegration
{
  public:
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          static IntegrationPoint point(
              Vector3(0.333333333333333,0.333333333333333,0),1.0/2.0);
          return &point;
        }
        virtual int getAccuracy() const {return 1;}
    };
    class N2 : public Integration
    {
      public:
        virtual int countPoints() const {return 3;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[3]=
{ IntegrationPoint(Vector3(0.666666666666667,0.166666666666667,0),0.333333333333333/2.0),
  IntegrationPoint(Vector3(0.166666666666667,0.666666666666667,0),0.333333333333333/2.0),
  IntegrationPoint(Vector3(0.166666666666667,0.166666666666667,0),0.333333333333333/2.0), };
          return points+i;
        }
        virtual int getAccuracy() const {return 2;}
    };
    class N3 : public Integration
    {
      public:
        virtual int countPoints() const {return 4;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[4]=
{ IntegrationPoint(Vector3(0.333333333333333,0.333333333333333,0),-0.562500000000000/2.0),
  IntegrationPoint(Vector3(0.600000000000000,0.200000000000000,0), 0.520833333333333/2.0),
  IntegrationPoint(Vector3(0.200000000000000,0.600000000000000,0), 0.520833333333333/2.0),
  IntegrationPoint(Vector3(0.200000000000000,0.200000000000000,0), 0.520833333333333/2.0) };
          return points+i;
        }
        virtual int getAccuracy() const {return 3;}
    };
    class N4 : public Integration
    {
      public:
        virtual int countPoints() const {return 6;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[6]=
{ IntegrationPoint(Vector3(0.816847572980459,0.091576213509771,0),0.109951743655322/2.0),
  IntegrationPoint(Vector3(0.091576213509771,0.816847572980459,0),0.109951743655322/2.0),
  IntegrationPoint(Vector3(0.091576213509771,0.091576213509771,0),0.109951743655322/2.0),
  IntegrationPoint(Vector3(0.108103018168070,0.445948490915965,0),0.223381589678011/2.0),
  IntegrationPoint(Vector3(0.445948490915965,0.108103018168070,0),0.223381589678011/2.0),
  IntegrationPoint(Vector3(0.445948490915965,0.445948490915965,0),0.223381589678011/2.0) };
          return points+i;
        }
        virtual int getAccuracy() const {return 4;}
    };
    class N5 : public Integration
    {
      public:
        virtual int countPoints() const {return 7;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[7]=
{ IntegrationPoint(Vector3(0.333333333333333,0.333333333333333,0),0.225000000000000/2.0),
  IntegrationPoint(Vector3(0.797426985353087,0.101286507323456,0),0.125939180544827/2.0),
  IntegrationPoint(Vector3(0.101286507323456,0.797426985353087,0),0.125939180544827/2.0),
  IntegrationPoint(Vector3(0.101286507323456,0.101286507323456,0),0.125939180544827/2.0),
  IntegrationPoint(Vector3(0.059715871789770,0.470142064105115,0),0.132394152788506/2.0),
  IntegrationPoint(Vector3(0.470142064105115,0.059715871789770,0),0.132394152788506/2.0),
  IntegrationPoint(Vector3(0.470142064105115,0.470142064105115,0),0.132394152788506/2.0) };
          return points+i;
        }
        virtual int getAccuracy() const {return 5;}
    };
    class N6 : public Integration
    {
      public:
        virtual int countPoints() const {return 12;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[12]=
{ IntegrationPoint(Vector3(0.873821971016996,0.063089014491502,0),0.050844906370207/2.0),
  IntegrationPoint(Vector3(0.063089014491502,0.873821971016996,0),0.050844906370207/2.0),
  IntegrationPoint(Vector3(0.063089014491502,0.063089014491502,0),0.050844906370207/2.0),
  IntegrationPoint(Vector3(0.501426509658179,0.249286745170910,0),0.116786275726379/2.0),
  IntegrationPoint(Vector3(0.249286745170910,0.501426509658179,0),0.116786275726379/2.0),
  IntegrationPoint(Vector3(0.249286745170910,0.249286745170910,0),0.116786275726379/2.0),
  IntegrationPoint(Vector3(0.636502499121399,0.310352451033785,0),0.082851075618374/2.0),
  IntegrationPoint(Vector3(0.636502499121399,0.053145049844816,0),0.082851075618374/2.0),
  IntegrationPoint(Vector3(0.310352451033785,0.636502499121399,0),0.082851075618374/2.0),
  IntegrationPoint(Vector3(0.310352451033785,0.053145049844816,0),0.082851075618374/2.0),
  IntegrationPoint(Vector3(0.053145049844816,0.636502499121399,0),0.082851075618374/2.0),
  IntegrationPoint(Vector3(0.053145049844816,0.310352451033785,0),0.082851075618374/2.0) };
          return points+i;
        }
        virtual int getAccuracy() const {return 6;}
    };
    virtual int countIntegrations() const {return 6;}
    virtual Integration const* getIntegration(int i) const
    {
      static N1 i1;
      static N2 i2;
      static N3 i3;
      static N4 i4;
      static N5 i5;
      static N6 i6;
      static Integration* integrations[6] = 
      {&i1,&i2,&i3,&i4,&i5,&i6};
      return integrations[i];
    }
};

class QuadIntegration : public EntityIntegration
{
  public:
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          static IntegrationPoint point(Vector3(0,0,0),4);
          return &point;
        }
        virtual int getAccuracy() const {return 1;}
    };
    class N2 : public Integration
    {
      public:
        virtual int countPoints() const {return 4;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static double const a =  0.577350269189626;
          static IntegrationPoint points[4]=
          { IntegrationPoint(Vector3(-a,-a,0),1),
            IntegrationPoint(Vector3( a,-a,0),1),
            IntegrationPoint(Vector3( a, a,0),1),
            IntegrationPoint(Vector3(-a, a,0),1) };
          return points + i;
        }
        virtual int getAccuracy() const {return 3;}
    };
    class N3 : public Integration
    {
      public:
        virtual int countPoints() const {return 9;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static double const a = 0.774596669241483;
          static double const b = 0.888888888888889;
          static double const c = 0.555555555555556;
          static IntegrationPoint points[9]=
          { IntegrationPoint(Vector3(-a,-a,0), c * c),
            IntegrationPoint(Vector3( 0,-a,0), b * c),
            IntegrationPoint(Vector3( a,-a,0), c * c),
            IntegrationPoint(Vector3(-a, 0,0), c * b),
            IntegrationPoint(Vector3( 0, 0,0), b * b),
            IntegrationPoint(Vector3( a, 0,0), c * b),
            IntegrationPoint(Vector3(-a, a,0), c * c),
            IntegrationPoint(Vector3( 0, a,0), b * c),
            IntegrationPoint(Vector3( a, a,0), c * c) };
          return points + i;
        }
        virtual int getAccuracy() const {return 5;}
    };
    virtual int countIntegrations() const {return 3;}
    virtual Integration const* getIntegration(int i) const
    {
      static N1 i1;
      static N2 i2;
      static N3 i3;
      static Integration* integrations[3] = 
      {&i1,&i2,&i3};
      return integrations[i];
    }
};

class TetrahedronIntegration : public EntityIntegration
{
  public:
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          static IntegrationPoint point(
              Vector3(0.25,0.25,0.25),1.0/6.0);
          return &point;
        }
        virtual int getAccuracy() const {return 1;}
    };
    class N2 : public Integration
    {
      public:
        virtual int countPoints() const {return 4;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[4]=
{ IntegrationPoint(Vector3(0.138196601125011,0.138196601125011,0.138196601125011),0.25/6.0), 
  IntegrationPoint(Vector3(0.585410196624969,0.138196601125011,0.138196601125011),0.25/6.0),
  IntegrationPoint(Vector3(0.138196601125011,0.585410196624969,0.138196601125011),0.25/6.0),
  IntegrationPoint(Vector3(0.138196601125011,0.138196601125011,0.585410196624969),0.25/6.0) };
          return points+i;
        }
        virtual int getAccuracy() const {return 2;}
    };
    class N3 : public Integration
    {
      public:
        virtual int countPoints() const {return 5;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static IntegrationPoint points[5]=
{ IntegrationPoint(Vector3(1./4., 1./4., 1./4.),-2./15.),
  IntegrationPoint(Vector3(1./6., 1./6., 1./6.),3./40.),
  IntegrationPoint(Vector3(1./6., 1./6., 1./2.),3./40.),
  IntegrationPoint(Vector3(1./6., 1./2., 1./6.),3./40.),
  IntegrationPoint(Vector3(1./2., 1./6., 1./6.),3./40.) };
          return points+i;
        }
        virtual int getAccuracy() const {return 3;}
    };
    virtual int countIntegrations() const {return 3;}
    virtual Integration const* getIntegration(int i) const
    {
      static N1 i1;
      static N2 i2;
      static N3 i3;
      static Integration* integrations[3] = 
      {&i1,&i2,&i3};
      return integrations[i];
    }
};

class PrismIntegration : public EntityIntegration
{
  public: /* a cheap one-point prism integration rule, not
             really legitimately 1st-order accurate but
             should be all right for our current purposes
           \todo revise this with some proper mathematics */
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          static IntegrationPoint point(Vector3(1.0/3.0,1.0/3.0,0),1);
          return &point;
        }
        virtual int getAccuracy() const {return 1;} //sort of
    };
    virtual int countIntegrations() const {return 1;}
    virtual Integration const* getIntegration(int) const
    {
      static N1 i1;
      return &i1;
    }
};

class PyramidIntegration : public EntityIntegration
{
  public:
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          /* one-point rule from
colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch12.d/AFEM.Ch12.pdf */
          static IntegrationPoint point(
              Vector3(0,0,-1.0/2.0),
              128.0 / 27.0);
          return &point;
        }
        virtual int getAccuracy() const {return 1;} //sort of
    };
    virtual int countIntegrations() const {return 1;}
    virtual Integration const* getIntegration(int) const
    {
      static N1 i1;
      return &i1;
    }
};

class HexahedronIntegration : public EntityIntegration
{
  public:
    class N1 : public Integration
    {
      public:
        virtual int countPoints() const {return 1;}
        virtual IntegrationPoint const* getPoint(int) const
        {
          static IntegrationPoint point(Vector3(0,0,0),8);
          return &point;
        }
        virtual int getAccuracy() const {return 1;}
    };
    class N2 : public Integration
    {
      public:
        virtual int countPoints() const {return 8;}
        virtual IntegrationPoint const* getPoint(int i) const
        {
          static double const x = 0.577350269189626;
          static IntegrationPoint points[8]=
          { IntegrationPoint(Vector3( x, x, x),1),
            IntegrationPoint(Vector3(-x, x, x),1),
            IntegrationPoint(Vector3( x,-x, x),1),
            IntegrationPoint(Vector3(-x,-x, x),1),
            IntegrationPoint(Vector3( x, x,-x),1),
            IntegrationPoint(Vector3(-x, x,-x),1),
            IntegrationPoint(Vector3( x,-x,-x),1),
            IntegrationPoint(Vector3(-x,-x,-x),1) };
          return points + i;
        }
        virtual int getAccuracy() const {return 3;}
    };
    virtual int countIntegrations() const {return 2;}
    virtual Integration const* getIntegration(int i) const
    {
      static N1 i1;
      static N2 i2;
      static Integration* integrations[2] = 
      {&i1,&i2};
      return integrations[i];
    }
};

EntityIntegration const* getIntegration(int meshEntityType)
{
  static EdgeIntegration edge;
  static TriangleIntegration triangle;
  static QuadIntegration quad;
  static TetrahedronIntegration tet;
  static PrismIntegration prism;
  static PyramidIntegration pyramid;
  static HexahedronIntegration hex;
  EntityIntegration* integrations[Mesh::TYPES] =
  {NULL,      //vertex
   &edge,     //edge
   &triangle, //triangle
   &quad,     //quad
   &tet,      //tet
   &hex,      //hex
   &prism,    //prism
   &pyramid}; //pyramid
  return integrations[meshEntityType];
}

Integrator::Integrator(int o):
  order(o)
{
}

Integrator::~Integrator()
{
}

void Integrator::inElement(MeshElement*)
{
}

void Integrator::outElement()
{
}

void Integrator::parallelReduce()
{
}

void Integrator::process(Mesh* m)
{
  int d = m->getDimension();
  MeshEntity* entity;
  MeshIterator* elements = m->begin(d);
  while ((entity = m->iterate(elements)))
  {
    if ( ! m->isOwned(entity)) continue;
    MeshElement* e = createMeshElement(m,entity);
    this->process(e);
    destroyMeshElement(e);
  }
  m->end(elements);
  this->parallelReduce();
}

void Integrator::process(MeshElement* e)
{
  this->inElement(e);
  int np = countIntPoints(e,this->order);
  for (int p=0; p < np; ++p)
  {
    Vector3 point;
    getIntPoint(e,this->order,p,point);
    double w = getIntWeight(e,this->order,p);
    double dV = getDV(e,point);
    this->atPoint(point,w,dV);
  }
  this->outElement();
}

class Measurer : public Integrator
{
  public:
    Measurer(int order):Integrator(order),m(0) {}
    void atPoint(Vector3 const&, double w, double dV)
    {
      m += w*dV;
    }
    double m;
};

double measure(MeshElement* e)
{
  Measurer measurer(getOrder(e));
  measurer.process(e);
  return measurer.m;
}

}//namespace apf
