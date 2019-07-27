#include "crvFaceOptim.h"
#include "LBFGS.h"
#include "crv.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvMath.h"
#include <iostream>
#include "apfMatrix.h"

namespace crv{

int CrvFaceReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

void CrvFaceReshapeObjFunc :: getInitFaceN()
{
  apf::Vector3 intFaceX;
  int numFNodes = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  for (int j = 0; j < numFNodes; j++) {
      mesh->getPoint(face, j, intFaceX);
      ifn.push_back(intFaceX);
  }
}

void CrvFaceReshapeObjFunc :: getInitTetN()
{
  apf::Adjacent adjT;
  apf::Vector3 intTetX;
  mesh->getAdjacent(face, 3, adjT);
  int numTNodes = mesh->getShape()->countNodesOn(mesh->TET);
  for (std::size_t i = 0; i < adjT.getSize(); i++) {
    for (int j = 0; j < numTNodes; j++) {
      mesh->getPoint(adjT[i], j, intTetX);
      itn.push_back(intTetX);
    }
  }
}

std::vector<double> CrvFaceReshapeObjFunc :: getInitialGuess()
{
  return convertNodeVectorToX(ifn);
}

std::vector<double> CrvFaceReshapeObjFunc :: convertNodeVectorToX(std::vector<apf::Vector3> fn)
{
  std::vector<double> x0;

  for (std::size_t i = 0; i < fn.size(); i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(fn[i][j]);
  
  return x0;
}

std::vector<apf::Vector3> CrvFaceReshapeObjFunc :: convertXtoNodeVector(const std::vector<double> &x)
{
  std::vector<apf::Vector3> a;
  apf::Vector3 v;
  std::size_t num = x.size()/d;
  //int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (std::size_t i = 0; i < num; i++) {
    v = {x[d*i], x[d*i + 1], x[d*i + 2]};
    a.push_back(v);
  }
  
  return a;
}

void CrvFaceReshapeObjFunc :: updateNodes(std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te)
{
  int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(face));
  for (int j = 0; j < numFNodes; j++)
    mesh->setPoint(face, j, fa[j]);

  if (d == 3 && P > 3) {
    apf::Adjacent adjT;
    mesh->getAdjacent(face, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      int numTNodes = mesh->getShape()->countNodesOn(mesh->getType(adjT[i]));
      for (int j =0; j < numTNodes; j++)
	mesh->setPoint(adjT[i], j, te[numTNodes*i + j]);
    }
  }
}

void CrvFaceReshapeObjFunc :: setNodes(std::vector<double> &x)
{
  std::vector<apf::Vector3> fn;// (ifn.begin(), ifn.end());
  std::vector<apf::Vector3> tn;// (itn.begin(), itn.end());
  std::vector<apf::Vector3> nod = convertXtoNodeVector(x);
  //blendTris(en, fn);
  //from all internal node vector to distict vector of nodes
  
  int nFN = mesh->getShape()->countNodesOn(mesh->getType(face));
  for (int i = 0; i <nFN; i++)
    fn.push_back(nod[i]);

  if (d > 2 && P > 3) {
    //int nTN = mesh->getShape()->countNodesOn(mesh->TET);
    for (std::size_t i = nFN; i <nod.size(); i++)
      tn.push_back(nod[i]);
  }

  updateNodes(fn, tn);
}

std::vector<double> CrvFaceReshapeObjFunc :: getVolume()
{
  if (d == 3) {
    apf::Adjacent adjT;
    apf::Matrix3x3 m;
    mesh->getAdjacent(face, 3, adjT);
    apf::Vector3 point0, point;
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      apf::Adjacent adjV;
      mesh->getAdjacent(adjT[i], 0, adjV);
      for (std::size_t j = 0; j < adjV.getSize(); j++) {
	if ( j == 0)
	  mesh->getPoint(adjV[j], 0, point0);
	else {
	  mesh->getPoint(adjV[j], 0, point);
	  for (int k = 0; k < 3; k++)
	    m[j-1][k] = point[k] - point0[k];
	}
      }
      double v = getDeterminant(m)/6.0;
      vol.push_back(v);
    }
  }
  return vol;
}

static double getAr(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2)
{
    p1 = p1 - p0;
    p2 = p2 - p0;
    double area = (apf::cross(p1, p2)).getLength();
    return (area/2.0);
}

double CrvFaceReshapeObjFunc :: computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm)
{
  int weight = 1;
  double sumf = 0;
  if (d == 3) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  for (int L = 0; L <= d*(P-1); L++) {
	    if ((I == J && J == K && I == 0) || (J == K && K == L && J == 0) || (I == K && K == L && I == 0) || (I == J && J == L && I == 0))
	      weight = 4;
	    else if ((I == J && I == 0) || (I == K && I == 0) || (I == L && I == 0) || (J == K && J == 0) || (J == L && J == 0) || (K == L && K == 0))
	      weight = 2;
	    else
	      weight = 1;
	    if (I + J + K + L == d*(P-1)) {
	      double fun = Nijkl(nodes,P,I,J,K)/(6.0*volm) - 1.0;
	      //std::cout<<"["<<I<<","<<J<<","<<K<<","<<L<<"]   "<<f<<std::endl;
	      sumf = sumf + weight*fun*fun;
	    }
	  }
	}
      }
    }
  }

  return sumf;
}

void CrvFaceReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ifn, itn);
}

double CrvFaceReshapeObjFunc :: getValue(std::vector<double> &x)
{
  setNodes(x);

  double sum = 0.0;

  if (d == 3) {
    apf::Adjacent adjT;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(face, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      apf::Element* el = apf::createElement(mesh->getCoordinateField(), adjT[i]);
      apf::getVectorNodes(el, allNodes);
      sum = sum + computeFValOfElement(allNodes, vol[i]);
      apf::destroyElement(el);
    }
/*
    apf::NewArray<apf::Vector3> fCP;
    apf::Vector3 xif;
    double a[3] = {1.0, 1.0, 1.0};
    double b[3] = {1.0, 1.0, 1.0};
    double aratio = 0.0;
    //double wfactor = 1.0;
    double gamma = 0.0;

    apf::Element* Fal = apf::createElement(mesh->getCoordinateField(), face);
    apf::getVectorNodes(Fal, fCP);
    int nFN = mesh->getShape()->countNodesOn(mesh->getType(face));

    int vN[3] = {getTriNodeIndex(P, P, 0), getTriNodeIndex(P, 0, P), getTriNodeIndex(P, 0, 0)};

    double triAphys = getAr(fCP[0], fCP[1], fCP[2]);
    apf::Vector3 prt0 = {1, 0, 0};
    apf::Vector3 prt1 = {0, 1, 0};
    apf::Vector3 prt2 = {0, 0, 1};
    double triAparnt = getAr(prt0, prt1, prt2);

    for (int j = 0; j < nFN; j++) {
      getBezierNodeXi(mesh->getType(face), P, j, xif);
      apf::Vector3 xifm = {1.0-xif[0]-xif[1], xif[0], xif[1]};
      a[0] = getAr(fCP[3*P + j], fCP[vN[0]], fCP[vN[1]]);
      b[0] = getAr(xifm, prt0, prt1);
      a[1] = getAr(fCP[3*P + j], fCP[vN[1]], fCP[vN[2]]);
      b[1] = getAr(xifm, prt1, prt2);
      a[2] = getAr(fCP[3*P + j], fCP[vN[2]], fCP[vN[0]]);
      b[2] = getAr(xifm, prt2, prt1);

      for (int k = 0; k < 3; k++) {
      	aratio = (a[k]*triAparnt/(b[k]*triAphys) - 1.0);
      	gamma = gamma + aratio*aratio;
      }
    }

    apf::destroyElement(Fal);

    sum = sum*(1 + gamma);
    */
    restoreInitialNodes();
  }
  return sum;
}

std::vector<double> CrvFaceReshapeObjFunc :: getGrad(std::vector<double> &x)
{
  //double fold = getValue(x);
  double eps = 1.0e-4;
  double h = eps;
  std::vector<double> g;
  double xmx = x[0];
  double xmn = x[0];
  double ymx = x[1];
  double ymn = x[1];
  double zmx = x[2];
  double zmn = x[2];
  double df = 0.0, dff = 0.0, dbb = 0.0;

  for (std::size_t i = 0; i < x.size(); i+=3) {
    if (x[i] >= xmx) xmx = x[i];
    if (x[i] <= xmn) xmn = x[i];
  }

  for (std::size_t i = 1; i < x.size(); i+=3) {
    if (x[i] >= ymx) ymx = x[i];
    if (x[i] <= ymn) ymn = x[i];
  }

  for (std::size_t i = 2; i < x.size(); i+=3) {
    if (x[i] >= zmx) zmx = x[i];
    if (x[i] <= zmn) zmn = x[i];
  }

  double delx = std::abs(xmx - xmn);
  double dely = std::abs(ymx - ymn);
  double delz = std::abs(zmx - zmn);
  double delta = 1.0;

  for (std::size_t i = 0; i < x.size(); i++) {
    if (i % 3 == 0) delta = delx;
    if (i % 3 == 1) delta = dely;
    if (i % 3 == 2) delta = delz;
    
    if (delta < eps)
      h = eps * std::abs(x[i]);
    else
      h = eps * delta;

    //if (std::abs(x[i]) > eps)
    //  h = eps * std::abs(x[i]);
    //else
    //  h = eps;

    x[i] = x[i] + h;
    double ff = getValue(x);
    x[i] = x[i] - h;
    
    x[i] = x[i] - h;
    double fb = getValue(x);
    x[i] = x[i] + h;
    
    df = (ff - fb)/(2.0 * h);
    g.push_back(df);
  }
  return g;
}

void CrvFaceOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvFaceOptim :: setTol(double t)
{
  tol = t;
}

bool CrvFaceOptim :: run()
{
  CrvFaceReshapeObjFunc *objF = new CrvFaceReshapeObjFunc(mesh, face);
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  LBFGS *l = new LBFGS(tol, iter, x0, objF);

  if (l->run()) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);
 /*
    apf::Adjacent adjT;
    mesh->getAdjacent(face, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      if (checkValidity(mesh, adjT[i], 1) > 1) {
      	objF->restoreInitialNodes();
	std::cout<<"invalid entity after fop with code--- "<<checkValidity(mesh, adjT[i], 1) <<std::endl;
      	//return false;
      }
    }
*/
    return true;
  }
  else {
    //std::cout<<"*****FaceOptim FAILURE"<<std::endl;
    return false;
  }
}

}
