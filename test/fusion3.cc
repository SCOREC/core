#include <ma.h>
#include <apf.h>
#include <gmi_analytic.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>
#include <vector>
using std::vector;
class Expression
{
public:
  virtual double eval(double x) const = 0;
};
class BSpline: public Expression
{
public:
  BSpline(int order_p, vector<double> &ctrlPts_p, vector<double> & knots_p);
  BSpline():order(-1){}
  ~BSpline () {};
  virtual double eval(double x) const;
private:
  int order;
  vector <double> ctrlPts;
  vector <double> knots;
};
BSpline :: BSpline(int order_p, vector<double> &ctrlPts_p, vector<double> & knots_p)
{
  assert(order_p>1);
  order=order_p;
  ctrlPts=ctrlPts_p;
  knots=knots_p;
}
double BSpline :: eval(double x) const
{
  // first find the interval of x in knots
  int leftKnot=order-1;
  int leftPt=0;
  while(knots.at(leftKnot+1)<x)
  {
    leftKnot++;
    leftPt++;
    if(leftKnot == ((int)knots.size()) - 1)
      break;
  }
  vector<double> pts(&(ctrlPts[leftPt]),&ctrlPts[leftPt+order]);
  vector<double> localKnots(&(knots[leftKnot-order+2]),&(knots[leftKnot+order]));
  for( int r=1; r<=order; r++)
  {
    // from bottom to top to save a buff
    for( int i=order-1; i>=r; i-- )
    {
      double a_left=localKnots.at(i-1);
      double a_right=localKnots.at(i+order-r-1);
      double alpha;
      if(a_right==a_left) alpha=0.; // not sure??
      else  alpha=(x-a_left)/(a_right-a_left);
      pts.at(i)=(1.-alpha)*pts.at(i-1)+alpha*pts.at(i);
    }
  }
  return pts.at(order-1);
}
void evalCoord(double para, double *xyz, void* userdata)
{
  Expression* R_p=((Expression**)userdata)[0];
  Expression* Z_p=((Expression**)userdata)[1];
  xyz[0]=R_p->eval(para);
  xyz[1]=Z_p->eval(para);
  //xyz[2]=0.0;
}
const int numGV=6;
const int numGE=6;
const int numGF=3;
double gvtx[6][3]={{1.016280, -0.491436, 0.000000},
                    {1.016036, 0.820123, 0.000000},
                    {0.996280, -0.491450, 0.000000},
                    {0.996036, 0.820093, 0.000000},
                    {0.500000, 0.000000, 0.000000},
                    {3.000000, 0.000000, 0.000000}};
double modelLen=2.5;
double center[]={1.7, 0.08,0.};
int splineOrder[6]={4,6,4,6,6,6};
double knots1[]={0.000000,0.000000,0.000000,0.000000,0.018182,0.036364,0.054545,0.072727,0.090909,0.109091,0.127273,0.145455,0.163636,0.181818,0.200000,0.218182,0.236364,0.254545,0.272727,0.290909,0.309091,0.327273,0.345455,0.363636,0.381818,0.400000,0.418182,0.436364,0.454545,0.472727,0.490909,0.509091,0.527273,0.545455,0.563636,0.581818,0.600000,0.618182,0.636364,0.654545,0.672727,0.690909,0.709091,0.727273,0.745455,0.763636,0.781818,0.800000,0.818182,0.836364,0.854545,0.872727,0.890909,0.909091,0.927273,0.945455,0.963636,0.981818,1.000000,1.000000,1.000000,1.000000};
double ctrPts1[]={1.016280,-0.491436,1.016339,-0.576055,1.016458,-0.745292,1.013586,-0.885286,1.005280,-1.187189,1.018980,-1.259631,1.145615,-1.358098,1.279885,-1.374623,1.371504,-1.366860,1.372193,-1.324668,1.374535,-1.254708,1.466892,-1.244305,1.579787,-1.249065,1.688906,-1.253107,1.765520,-1.250035,1.767502,-1.216165,1.786996,-1.173692,1.894360,-1.111232,2.019187,-1.050151,2.133636,-0.971656,2.214745,-0.791719,2.292841,-0.590573,2.361567,-0.392790,2.370055,-0.221740,2.366193,-0.125585,2.364867,-0.042481,2.364982,0.042260,2.365556,0.127201,2.367961,0.226271,2.365545,0.384391,2.327108,0.519472,2.248071,0.724465,2.147787,0.983202,2.049837,1.059713,1.796218,1.079554,1.643892,1.075869,1.596925,1.094105,1.529310,1.157013,1.445690,1.228385,1.374921,1.292198,1.371379,1.325309,1.362306,1.348841,1.288830,1.350508,1.275925,1.343253,1.278241,1.330827,1.264274,1.290961,1.242774,1.244102,1.222321,1.212539,1.202955,1.194695,1.180196,1.181590,1.137059,1.168767,1.094365,1.164451,1.044741,1.164281,1.014930,1.163964,1.013266,1.120499,1.015576,1.122973,1.015882,0.921073,1.016036,0.820123};
double knots2[]={0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
double ctrPts2[]={1.016036,0.820123,1.016096,0.780123,1.016157,0.740123,1.016252,-0.451436,1.016252,-0.451436,1.016280,-0.491436};
double knots3[]={0.000000,0.000000,0.000000,0.000000,0.018182,0.036364,0.054545,0.072727,0.090909,0.109091,0.127273,0.145455,0.163636,0.181818,0.200000,0.218182,0.236364,0.254545,0.272727,0.290909,0.309091,0.327273,0.345455,0.363636,0.381818,0.400000,0.418182,0.436364,0.454545,0.472727,0.490909,0.509091,0.527273,0.545455,0.563636,0.581818,0.600000,0.618182,0.636364,0.654545,0.672727,0.690909,0.709091,0.727273,0.745455,0.763636,0.781818,0.800000,0.818182,0.836364,0.854545,0.872727,0.890909,0.909091,0.927273,0.945455,0.963636,0.981818,1.000000,1.000000,1.000000,1.000000};
double ctrPts3[]={0.996280,-0.491450,0.996314,-0.575982,0.996384,-0.745045,0.993883,-0.885439,0.984202,-1.183789,1.003006,-1.274805,1.137848,-1.376951,1.278533,-1.393845,1.380030,-1.391029,1.396554,-1.314309,1.388522,-1.275221,1.464230,-1.264109,1.579747,-1.269291,1.686975,-1.272306,1.772659,-1.273010,1.791923,-1.213667,1.797690,-1.193204,1.903648,-1.128151,2.027706,-1.068902,2.150735,-0.983387,2.233518,-0.798453,2.311373,-0.598166,2.381278,-0.397133,2.390131,-0.221365,2.386160,-0.124820,2.384875,-0.042445,2.384981,0.042215,2.385553,0.126855,2.387960,0.225757,2.385552,0.386796,2.345935,0.526927,2.266178,0.731424,2.168455,0.991192,2.052354,1.081869,1.797203,1.098672,1.642207,1.097144,1.611419,1.109567,1.542371,1.171822,1.458626,1.243422,1.389165,1.307523,1.396699,1.321986,1.363915,1.372893,1.292112,1.372394,1.253441,1.351408,1.259218,1.333175,1.245723,1.299509,1.224933,1.252985,1.206429,1.225063,1.190934,1.211020,1.173123,1.200473,1.133433,1.188578,1.092442,1.183999,1.050234,1.186137,0.994144,1.176988,0.993480,1.116593,0.995524,1.123680,0.995865,0.921289,0.996036,0.820093};
double knots4[]={0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
double ctrPts4[]={0.996036,0.820093,0.996096,0.780093,0.996157,0.740093,0.996252,-0.451450,0.996252,-0.451450,0.996280,-0.491450};
double knots5[]={0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
double ctrPts5[]={0.500000,0.000000,0.500000,-2.200000,0.500000,-3.400000,3.000000,-2.400000,3.000000,-1.200000,3.000000,0.000000};
double knots6[]={0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
double ctrPts6[]={3.000000,0.000000,3.000000,1.200000,3.000000,2.400000,0.500000,3.400000,0.500000,2.200000,0.500000,0.000000};
void makeBSpline(BSpline** splines, int  order, int numPts, double* ctrlPts, double * knots)
{
  vector<double> X_p, Y_p, knots_p;
  for (int i=0; i<numPts;i++)
  {
    X_p.push_back(ctrlPts[2*i]);
    Y_p.push_back(ctrlPts[2*i+1]);
  }
  for( int i=0; i<numPts+order; i++)
  {
    knots_p.push_back(knots[i]);
  }

  splines[0] = new BSpline(order,X_p,knots_p);
  splines[1] = new BSpline(order,Y_p,knots_p);
}

void edgeFunction(double const p[2], double *xyz, void*  data)
{
  evalCoord(p[0], xyz, data);
  xyz[2]=0.;
}

void faceFunction(double const p[2], double x[3], void * data)
{
  (void)p;(void)x;(void)data;
}

void vertexFunction(double const p[2], double x[3], void * data)
{
  (void)p;(void)x;(void)data;
}

void reparam_zero(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 0;
  to[1] = 0;
}

void reparam_one(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 1;
  to[1] = 0;
}

agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

void make_edge_topo(gmi_model* m, gmi_ent* e, int v0tag, int v1tag)
{
  agm_bdry b = add_bdry(m, e);
  agm_use u0 = add_adj(m, b, v0tag);
  gmi_add_analytic_reparam(m, u0, reparam_zero, 0);
  agm_use u1 = add_adj(m, b, v1tag);
  gmi_add_analytic_reparam(m, u1, reparam_one, 0);
}

int const edge_verts[numGE][2] = {
  {1,2},
  {2,1},
  {3,4},
  {4,3},
  {5,6},
  {6,5},
};

void make_faces_topo(gmi_model* m, gmi_ent* faces[numGF])
{
  agm* topo = gmi_analytic_topo(m);
  agm_bdry b = add_bdry(m, faces[0]);
  add_adj(m, b, 1);
  add_adj(m, b, 2);
  printf("face 1 has %d\n",agm_down_count(topo, agm_from_gmi(faces[0])));
  b = add_bdry(m, faces[1]);
  add_adj(m, b, 1);
  add_adj(m, b, 2);
  b = add_bdry(m, faces[1]);
  add_adj(m, b, 3);
  add_adj(m, b, 4);
  b = add_bdry(m, faces[2]);
  add_adj(m, b, 3);
  add_adj(m, b, 4);
  b = add_bdry(m, faces[2]);
  add_adj(m, b, 5);
  add_adj(m, b, 6);
}

gmi_model* makeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edgePeriodic = 0;
  double edgeRange[2] = {0, 1};
  for(int i=0; i<numGV; i++)
    gmi_add_analytic(model, 0, i+1, vertexFunction, NULL, NULL, NULL);
  double *knotsVec[]={knots1, knots2, knots3, knots4, knots5, knots6};
  double *ctrPtsVec[]={ctrPts1, ctrPts2,ctrPts3,ctrPts4,ctrPts5,ctrPts6};
  int numPts[]={sizeof(ctrPts1), sizeof(ctrPts2),sizeof(ctrPts3),sizeof(ctrPts4),sizeof(ctrPts5),sizeof(ctrPts6)};
  for(int i=0; i<numGV; i++)
    numPts[i]/=sizeof(double)*2;
  for(int i=0; i<numGE; i++)
  {
    BSpline** data= new BSpline*[2];
    makeBSpline(data, splineOrder[i], numPts[i], ctrPtsVec[i], knotsVec[i]);
    gmi_ent* edge =
      gmi_add_analytic(model, 1, i+1, edgeFunction, &edgePeriodic, &edgeRange, data);
    make_edge_topo(model, edge, edge_verts[i][0], edge_verts[i][1]);
  }
  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,0},{0,0}};
  gmi_ent* faces[numGF];
  for(int i=0; i<numGF; i++)
    faces[i] =
      gmi_add_analytic(model, 2, i+1, faceFunction, facePeriodic, faceRanges, 0);
  make_faces_topo(model, faces);
  return model;
}
class Vortex : public ma::AnisotropicFunction
{
  public:
    Vortex(ma::Mesh* m, double center[3], double len)
    {
      std::cout<<" Vortex AnisotropicFunction center len "<<center[0]<<" "<<center[1]<<" "<<center[2]<<" "<<len<<std::endl;
      mesh = m;
      modelLen = len;
      average = ma::getAverageEdgeLength(m);
      ma::Vector lower,upper;
      for(int i=0; i<3; i++)
        centroid[i]=center[i];
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
        dx=1.;
        dy=0;
      }
      double small=0.1*average;
      double large=0.3*average;
      h[0]=small+large*modelLen*fabs(r-modelLen/1.8);
      h[1]=large+large*fabs(r-modelLen/1.8);
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
    double modelLen;
    ma::Vector centroid;
};

int main(int argc, char * argv[])
{
  assert(argc==2);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_model* model = makeModel();
  gmi_write_dmg(model, "made.dmg");
  apf::Mesh2* mesh=apf::loadMdsMesh(model, argv[1]);
  mesh->verify();
  Vortex sfv(mesh, center, modelLen);
  ma::Input* in = ma::configure(mesh,&sfv);
  in->maximumIterations = 9;
  ma::adapt(in);
  mesh->verify();
  apf::writeVtkFiles("adapted",mesh);
  //clean data
  // to do
}
