#include <ma.h>
#include <apf.h>
#include <gmi_analytic.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <parma.h>
#include <apfZoltan.h>
#include <pcu_util.h>

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

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

gmi_model* makeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edgePeriodic = 1;
  double edgeRange[2] = {0, 2 * apf::pi};
  gmi_add_analytic(model, 1, 1, edgeFunction, &edgePeriodic, &edgeRange, 0);
  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,0},{0,0}};
  gmi_ent* face = gmi_add_analytic(model, 2, 1, faceFunction, facePeriodic, faceRanges, 0);


  agm_bdry b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(face));
  agm_use faceUse = add_adj(model, b, 1);
  gmi_add_analytic_reparam(model, faceUse, faceFunction, 0);

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
      PCU_ALWAYS_ASSERT( apf::getMdsIndex(m, e) == i );
      PCU_ALWAYS_ASSERT( apf::getMdsEntity(m, d, i) == e );
      ++i;
    }
    m->end(it);
  }
}

static void fusionAdapt(apf::Mesh2* m)
{
  Vortex sf(m);
  const ma::Input* in = ma::configure(m, &sf);
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
      mesh = apf::loadMdsMesh(model, meshFile, PCUObj.get());
      mesh->verify();
      testIndexing(mesh);
      fusionAdapt(mesh);
    } else {
      mesh = apf::makeEmptyMdsMesh(model, 2, false, PCUObj.get());
    }
  }
};

int main( int argc, char* argv[])
{
  PCU_ALWAYS_ASSERT(argc==2);
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  GroupCode code;
  code.model = makeModel();
  code.meshFile = argv[1];
  apf::Unmodulo outMap(pcu_obj.Self(), 2);
  Parma_SplitPartition(nullptr, 2, code, &pcu_obj);
  //Have to call switchPCU here because the mesh needed to be made inside GroupCode run()
  //and inside parma_group.cc runInGroups() we set code.PCUObj to groupedPCU 
  code.mesh->switchPCU(&pcu_obj);
  apf::remapPartition(code.mesh, outMap);
  code.mesh->verify();
  code.mesh->destroyNative();
  apf::destroyMesh(code.mesh);
  }
  MPI_Finalize();
}
