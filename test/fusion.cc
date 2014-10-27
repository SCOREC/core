#include <ma.h>
#include <apf.h>
#include <gmi_analytic.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>

double const a_param = 0.2;
double const b_param = 1.0;
double const c_param = 0.2;
double const d_param = 0.0;
double const e_param = 1.3;

void edgeFunction(double const p[2], double x[3], void*)
{
  double phi = p[0];
  x[0] = a_param + b_param*(cos(phi + c_param*sin(phi)));
  x[1] = d_param + e_param*sin(phi);
  x[2] = 0;
}

void faceFunction(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}

gmi_model* makeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edgePeriodic = 1;
  double edgeRange[2] = {0, 2 * apf::pi};
  gmi_add_analytic(model, 1, 1, edgeFunction, &edgePeriodic, &edgeRange, 0);
  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,0},{0,0}};
  gmi_add_analytic(model, 2, 1, faceFunction, facePeriodic, faceRanges, 0);
  return model;
}

class Vortex : public ma::AnisotropicFunction
{
  public:
    Vortex(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      ma::Vector lower,upper;
      lower[0]=a_param - b_param;
      lower[1]=d_param - e_param;
      lower[2]=0;
      upper[0]=a_param + b_param;
      upper[1]=d_param + e_param;
      upper[2]=0;
      for(int i=0; i<3; i++)
        centroid[i]=0.5*(upper[i]+lower[i]);
    }
    virtual void getValue(
        ma::Entity* v,
        ma::Matrix& R,
        ma::Vector& h)
    {
      ma::Vector x = ma::getPosition(mesh,v);
      double dx=x[0]-centroid[0];
      double dy=x[1]-centroid[1];
      double r=sqrt(dx*dx+dy*dy);
      if(r>1e-6) // if the vertex near to origin
      {
        dx=dx/r;
        dy=dy/r;
      }
      else
      {
        dx=1.0;
        dy=0.0;
      }
      double modelLen=b_param;
      h[0]=0.01+0.1*fabs(r-modelLen/3.);
      h[1]=0.06+0.1*fabs(r-modelLen/3.);
      h[2]=1.;
      R[0][0]=dx;
      R[1][0]=dy;
      R[2][0]=0;
      R[0][1]=-1.*dy;
      R[1][1]=dx;
      R[2][1]=0;   
      R[0][2]=0;
      R[1][2]=0;
      R[2][2]=1.;
    }
  private:
    ma::Mesh* mesh;
    double average;
    ma::Vector centroid;
};

static void testIndexing(apf::Mesh2* m)
{
  for (int d = 0; d <= m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    int i = 0;
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      assert( apf::getMdsIndex(m, e) == i );
      assert( apf::getMdsEntity(m, d, i) == e );
      ++i;
    }
    m->end(it);
  }
}

static void fusionAdapt(apf::Mesh2* m)
{
  Vortex sf(m);
  ma::Input* in = ma::configure(m, &sf);
  in->maximumIterations = 9;
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  ma::adapt(in);
  m->verify();
}

struct GroupCode : public Parma_GroupCode
{
  apf::Mesh2* mesh;
  gmi_model* model;
  const char* meshFile;
  void run(int group)
  {
    if (group == 0) {
      mesh = apf::loadMdsMesh(model, meshFile);
      mesh->verify();
      testIndexing(mesh);
      fusionAdapt(mesh);
    } else {
      mesh = apf::makeEmptyMdsMesh(model, 2, false);
    }
  }
};

int main( int argc, char* argv[])
{
  assert(argc==2);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  GroupCode code;
  code.model = makeModel();
  code.meshFile = argv[1];
  apf::Unmodulo outMap(PCU_Comm_Self(), 2);
  Parma_SplitPartition(NULL, 2, code);
  apf::remapPartition(code.mesh, outMap);
  code.mesh->verify();
  code.mesh->destroyNative();
  apf::destroyMesh(code.mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
