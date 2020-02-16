#include <apf.h>
#include "crv.h"
#include "crvSnap.h"
#include "crvObjectiveFunctions.h"
#include "crvBezier.h"
#include "crvQuality.h"
#include "crvMath.h"
#include "apfIntegrate.h"



static double getAr(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2)
{
  p1 = p1 - p0;
  p2 = p2 - p0;
  double area = (apf::cross(p1, p2)).getLength();
  return (area/2.0);
}

/* static double getLinearVolPhys(apf::Mesh2* m, apf::MeshEntity* e) */
/* { */
/*   apf::MeshEntity* vs[12]; */
/*   int n = m->getDownward(e, 0, vs); */
/*   apf::Vector3 coords[12]; */
/*   for (int i = 0; i < n; i++) { */
/*     m->getPoint(vs[i], 0, coords[i]); */
/*   } */

/*   if (m->getType(e) == apf::Mesh::TRIANGLE) */
/*   { */
/*     return apf::cross(coords[1]-coords[0], coords[2]-coords[0]).getLength()/2.; */
/*   } */
/*   else if (m->getType(e) == apf::Mesh::TET) */
/*   { */
/*     apf::Matrix3x3 J; */
/*     J[0] = coords[1] - coords[0]; */
/*     J[1] = coords[2] - coords[0]; */
/*     J[2] = coords[3] - coords[0]; */
/*     return apf::getDeterminant(J) / 6.; */
/*   } */
/*   else */
/*     PCU_ALWAYS_ASSERT_VERBOSE(0, */
/*     	"Not implemented for entities of type other than tri or tet!"); */
/*   return 0.; */
/* } */

static double computeFValNIJKL(apf::Mesh2* m, apf::MeshEntity* e, ma::SizeField* s = 0)
{
  PCU_ALWAYS_ASSERT_VERBOSE(s == 0, "Not implemented for non-zero sizefield!");

  int d = m->getDimension();
  int P = m->getShape()->getOrder();

  apf::NewArray<apf::Vector3> nodes;
  apf::Element* el = apf::createElement(m->getCoordinateField(), e);
  apf::getVectorNodes(el, nodes);
  apf::destroyElement(el);

  double volm = getLinearVolPhys(m, e);

  int weight = 1;
  double sumf = 0;
  if (d == 3) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  for (int L = 0; L <= d*(P-1); L++) {
	    if ((I == J && J == K && I == 0) ||
	    	(J == K && K == L && J == 0) ||
	    	(I == K && K == L && I == 0) ||
	    	(I == J && J == L && I == 0))
	      weight = 4;
	    else if ((I == J && I == 0) ||
	    	     (I == K && I == 0) ||
	    	     (I == L && I == 0) ||
	    	     (J == K && J == 0) ||
	    	     (J == L && J == 0) ||
	    	     (K == L && K == 0))
	      weight = 2;
	    else
	      weight = 1;
	    if (I + J + K + L == d*(P-1)) {
	      double f = crv::Nijkl(nodes,P,I,J,K)/(6.0*volm) - 1.0;
	      sumf = sumf + weight*f*f;
	    }
	  }
	}
      }
    }
  }

  if (d == 2) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  if ((I == J && I == 0) ||
	      (J == K && J == 0) ||
	      (I == K && I == 0))
	    weight = 2;
	  else
	    weight = 1;
	  if (I + J + K == d*(P-1)) {
	    double f = crv::Nijk(nodes,P,I,J)/(4.0*volm) - 1.0;
	    sumf = sumf + weight*f*f;
	  }
	}
      }
    }
  }
  return sumf;
}

static double computeFValDetJ(apf::Mesh2* m, apf::MeshEntity* e, ma::SizeField* s)
{
  apf::EntityIntegration *eInt = new apf::EntityIntegration();
  apf::Integration *ir = eInt->getIntegration(m->getType(e));
  int nip = eInt->countIntegrations();

  apf::MeshElement* me = apf::createMeshElement(m, e);
  Matrix3x3 J, T, Jm;

  double jDet, sum = 0.;

  for (int i = 0; i < nip; i++) {
    IntegrationPoint *ip = ir->getPoint(i);
    Vector3 qp = ip->param;
    double w = ip->weight;
    apf::getJacobian(me, qp, J);

    s->getTransform(me, qp, T);
    Jm = J*T; // Jacobian in metric space
    determinant(Jm, jDet);
    sum += (jDet - 1.) * (jDet - 1);
  }
  return sum;
}


namespace crv
{

int InternalEdgeReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

double InternalEdgeReshapeObjFunc :: getValue(const vector<double> &x)
{
  setNodes(x);
  double sum = 0.0;
  if (d == 2) {
    apf::Adjacent adjF;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 2, adjF);
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      sum = sum + computeFValNIJKL(mesh, adjF[i]);
    }
    restoreInitialNodes();
  }
  if (d == 3) {
    apf::Adjacent adjT;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      sum = sum + computeFValNIJKL(mesh, adjT[i]);
    }
    // TODO: In the original code this line is commented. Why?
    restoreInitialNodes();
  }
  return sum;
}


/* vector<double> InternalEdgeReshapeObjFunc :: getGrad(const vector<double> &_x) */
/* { */
/*   vector<double> x; */
/*   for (int i = 0; i < _x.size(); i++) */
/*     x.push_back(_x[i]); */

/*   // TODO: Why these values? */
/*   double eps = 1.0e-5; */
/*   double delta = 1.0; */
/*   double h = eps; */

/*   vector<double> g; */
/*   double xmx = x[0]; */
/*   double xmn = x[0]; */
/*   double ymx = x[1]; */
/*   double ymn = x[1]; */
/*   double zmx = x[2]; */
/*   double zmn = x[2]; */
/*   double df = 0.0, dff = 0.0, dbb = 0.0; */

/*   for (std::size_t i = 0; i < x.size(); i+=3) { */
/*     if (x[i] >= xmx) xmx = x[i]; */
/*     if (x[i] <= xmn) xmn = x[i]; */
/*   } */

/*   for (std::size_t i = 1; i < x.size(); i+=3) { */
/*     if (x[i] >= ymx) ymx = x[i]; */
/*     if (x[i] <= ymn) ymn = x[i]; */
/*   } */

/*   for (std::size_t i = 2; i < x.size(); i+=3) { */
/*     if (x[i] >= zmx) zmx = x[i]; */
/*     if (x[i] <= zmn) zmn = x[i]; */
/*   } */

/*   double delx = std::abs(xmx - xmn); */
/*   double dely = std::abs(ymx - ymn); */
/*   double delz = std::abs(zmx - zmn); */

/*   for (std::size_t i = 0; i < x.size(); i++) { */
/*     if (i % 3 == 0) delta = delx; */
/*     if (i % 3 == 1) delta = dely; */
/*     if (i % 3 == 2) delta = delz; */

/*     h = eps * delta; */

/*     x[i] = x[i] + h; */
/*     double ff = getValue(x); */
/*     x[i] = x[i] - h; */

/*     x[i] = x[i] - h; */
/*     double fb = getValue(x); */
/*     x[i] = x[i] + h; */

/*     df = (ff - fb)/(2.0 * h); */
/*     g.push_back(df); */
/*   } */
/*   return g; */
/* } */

vector<double> InternalEdgeReshapeObjFunc :: getInitialGuess()
{
  return convertNodeVectorToX(ien, ifn, itn);
}

void InternalEdgeReshapeObjFunc :: setNodes(const vector<double> &x)
{
  std::vector<apf::Vector3> en;// (ien.begin(), ien.end());
  std::vector<apf::Vector3> fn;// (ifn.begin(), ifn.end());
  std::vector<apf::Vector3> tn;// (itn.begin(), itn.end());
  std::vector<apf::Vector3> nod = convertXtoNodeVector(x);
  //blendTris(en, fn);
  //from all internal node vector to distict vector of nodes

  int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (int i = 0; i <nEN; i++)
    en.push_back(nod[i]);

  int nFN = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);

  for (std::size_t i = nEN; i < nEN + nFN*adjF.getSize(); i++)
    fn.push_back(nod[i]);

  if (d > 2 && P > 3) {
    //int nTN = mesh->getShape()->countNodesOn(mesh->TET);
    for (std::size_t i = nEN+nFN; i <nod.size(); i++)
      tn.push_back(nod[i]);
  }

  updateNodes(en, fn, tn);
}

void InternalEdgeReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ien, ifn, itn);
}

void InternalEdgeReshapeObjFunc :: getInitEdgeN()
{
  apf::Vector3 intEdgeX;
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));

  for (int i = 0; i < numENodes; i++) {
    mesh->getPoint(edge, i, intEdgeX);
    ien.push_back(intEdgeX);
  }
}

void InternalEdgeReshapeObjFunc :: getInitFaceN()
{
  apf::Adjacent adjF;
  apf::Vector3 intFaceX;
  mesh->getAdjacent(edge, 2, adjF);
  int numFNodes = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    for (int j = 0; j < numFNodes; j++) {
      mesh->getPoint(adjF[i], j, intFaceX);
      ifn.push_back(intFaceX);
    }
  }
}

void InternalEdgeReshapeObjFunc :: getInitTetN()
{
  apf::Adjacent adjT;
  apf::Vector3 intTetX;
  mesh->getAdjacent(edge, 3, adjT);
  int numTNodes = mesh->getShape()->countNodesOn(mesh->TET);
  for (std::size_t i = 0; i < adjT.getSize(); i++) {
    for (int j = 0; j < numTNodes; j++) {
      mesh->getPoint(adjT[i], j, intTetX);
      itn.push_back(intTetX);
    }
  }
}

vector<double> InternalEdgeReshapeObjFunc :: convertNodeVectorToX(
    const vector<apf::Vector3> &en,
    const vector<apf::Vector3> &fn,
    const vector<apf::Vector3> &tn)
{
  vector<double> x0;
  for (int i = 0; i < P-1; i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(en[i][j]);

  for (std::size_t i = 0; i < fn.size(); i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(fn[i][j]);

  // TODO: Why Is There a Need For if conditions?
  // Isn't tn.size() == 0 when there are no tets in the mesh?
  // or when the order of the mesh is less than 4
  if (d > 2 && P > 3) {
    for (std::size_t i = 0; i < tn.size(); i++)
      for (int j = 0; j < 3; j++)
	x0.push_back(tn[i][j]);
  }
  return x0;
}

vector<apf::Vector3> InternalEdgeReshapeObjFunc :: convertXtoNodeVector(const vector<double> &x)
{
  std::vector<apf::Vector3> a;
  apf::Vector3 v;
  if (d == 2) {
    std::size_t num = x.size()/d; // TODO: this seems unsafe
    //int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    //check later for 2D case:: x should not include z coordinate in optim search
    for (std::size_t i = 0; i < num; i++) {
      v = {x[d*i], x[d*i + 1], 0.0};
      a.push_back(v);
    }
  }
  if (d == 3) {
    std::size_t num = x.size()/d; // TODO: this seems unsafe
    //int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (std::size_t i = 0; i < num; i++) {
      v = {x[d*i], x[d*i + 1], x[d*i + 2]};
      a.push_back(v);
    }
  }
  return a;
}

void InternalEdgeReshapeObjFunc :: blendTris(
    const vector<apf::Vector3> &egn,
    vector<apf::Vector3> &faceNodes)
{
  apf::Vector3 xi;
  apf::Adjacent adjF;
  apf::Adjacent adjE;
  std::vector<apf::Vector3> cien (ien.begin(),ien.end());

  mesh->getAdjacent(edge, 2, adjF);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    mesh->getAdjacent(adjF[i], 1, adjE);
    int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
    for (std::size_t j = 0; j < adjE.getSize(); j++) {
      if (adjE[j] == edge) {
	int jj = 1;
	if ( j == 0)
	  jj = 2;
	else if ( j == 1)
	  jj = 0;
	else
	  jj = 1;

	for (int k = 0; k < numFNodes; k++) {
	  getBezierNodeXi(mesh->TRIANGLE, P, k, xi);
	  for (std::size_t ii = 0; ii < egn.size(); ii++) {
	    double factor =
	      0.5 * (xi[j]/(1-xi[jj])) * binomial(P, ii+1) * intpow(1-xi[jj], ii+1) * intpow(xi[jj], P-ii-1)
	    + 0.5 * (xi[jj]/(1-xi[j])) * binomial(P, ii+1) * intpow(xi[j], ii+1) * intpow(1-xi[j], P-ii-1);
	    faceNodes[numFNodes*i+k] = faceNodes[numFNodes*i+k] + (egn[ii] - cien[ii]) * factor;
	  }
	}
      }
    }
  }
}

void InternalEdgeReshapeObjFunc :: updateNodes(
    const vector<apf::Vector3> &ed,
    const vector<apf::Vector3> &fa,
    const vector<apf::Vector3> &te)
{
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (int i = 0; i < numENodes; i++)
    mesh->setPoint(edge, i, ed[i]);

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
    for (int j = 0; j < numFNodes; j++)
      mesh->setPoint(adjF[i], j, fa[numFNodes*i + j]);
  }

  if (d == 3 && P > 3) {
    apf::Adjacent adjT;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      int numTNodes = mesh->getShape()->countNodesOn(mesh->getType(adjT[i]));
      for (int j =0; j < numTNodes; j++)
  mesh->setPoint(adjT[i], j, te[numTNodes*i + j]);
    }
  }
}





// Boundary Edge Objective Functions Impl

int BoundaryEdgeReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

double BoundaryEdgeReshapeObjFunc :: getValue(const vector<double> &x)
{
  setNodes(x);

  double sum = 0.0;
  if (d == 2) {
    apf::Adjacent adjF;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 2, adjF);
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      sum = sum + computeFValNIJKL(mesh, adjF[i]);
    }
    restoreInitialNodes();
  }

  if (d == 3) {
    apf::Adjacent adjT;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      sum = sum + computeFValNIJKL(mesh, adjT[i]);
    }

    double ad = 0.0;
    double xr = 1.0;
    double xir = 1.0;
    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
    double beta = 0.0;
    std::vector<apf::Vector3> xs;
    xs.clear();
    std::vector<apf::Vector3> xis;
    xis.clear();
    xis.push_back(apf::Vector3(-1.0, 0.0, 0.0));
    for (int i = 0; i <nEN; i++) {
      apf::Vector3 currentXi;
      getBezierNodeXi(mesh->getType(edge), P, i, currentXi);
      xis.push_back(currentXi);
    }
    xis.push_back(apf::Vector3(+1.0, 0.0, 0.0));

    apf::MeshElement* me = apf::createMeshElement(mesh, edge);
    for (int i = 0; i < xis.size(); i++) {
      apf::Vector3 scord;
      apf::mapLocalToGlobal(me, xis[i], scord);
      xs.push_back(scord);
    }
    apf::destroyMeshElement(me);

    for (int i = 1; i < xs.size()-1; i++) {

      xr = (xs[i] - xs[0]).getLength() /
	   (xs[xs.size()-1] - xs[0]).getLength();
      xir = (xis[i] - xis[0]).getLength() /
	    (xis[xs.size()-1] - xis[0]).getLength();
      ad = (1.0*xr/xir - 1.0);   //(alpha*alpha);
      beta = beta + ad*ad;
      //sum = sum + ad*ad;
    }

    //sum = sum*(1 + beta);
    /* apf::destroyElement(Edl); */

    apf::Adjacent adjF;
    mesh->getAdjacent(edge, 2, adjF);

    std::vector<apf::Vector3> xfs;
    std::vector<apf::Vector3> xifs;
    double arPhys[3] = {1.0, 1.0, 1.0};
    double arParnt[3] = {1.0, 1.0, 1.0};
    double adf;
    double gamma = 0.0;

    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      xfs.clear();
      xifs.clear();
      adf = 0.0;

      if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
	int nFN = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
	apf::Vector3 faceXi;
	for (int j = 0; j < nFN; j++) {
	  getBezierNodeXi(mesh->getType(adjF[i]), P, j, faceXi);
	  xifs.push_back(faceXi);
	}
	xifs.push_back(apf::Vector3(0.0, 0.0, 0.0));
	xifs.push_back(apf::Vector3(1.0, 0.0, 0.0));
	xifs.push_back(apf::Vector3(0.0, 1.0, 0.0));

	apf::MeshElement* mef = apf::createMeshElement(mesh, adjF[i]);
	for (size_t k = 0; k < xifs.size(); k++) {
	  apf::Vector3 fcord;
	  apf::mapLocalToGlobal(mef, xifs[k], fcord);
	  xfs.push_back(fcord);
	}
	apf::destroyMeshElement(mef);

	double triPhys = getAr(xfs[xifs.size()-1], xfs[xifs.size()-2], xfs[xifs.size()-3]);
	double triParnt = getAr(xifs[xifs.size()-1], xifs[xifs.size()-2], xifs[xifs.size()-3]);

	for (int j = 0; j < nFN; j++) {
	  for (int k = 0; k < 2; k++) {
	    arPhys[k] = getAr(xfs[j], xfs[xifs.size()-(3-k)], xfs[xifs.size()-(2-k)]);
	    arParnt[k] = getAr(xifs[j], xifs[xifs.size()-(3-k)], xifs[xifs.size()-(2-k)]);
	    adf = (1.0*arPhys[k]*triParnt/(arParnt[k]*triPhys) - 1);
	    gamma = gamma + adf*adf;
	  }
	}
      }
    }
    /* double gamma = 0.0; */
    sum = sum*(1 + beta + 0.3*gamma);
    restoreInitialNodes();
  }
  return sum;
}

/* vector<double> BoundaryEdgeReshapeObjFunc :: getGrad(const vector<double> &_x) */
/* { */
/*   vector<double> x; */
/*   for (int i = 0; i < _x.size(); i++) { */
/*     x[i] = _x[i]; */
/*   } */

/*   //double fold = getValue(x); */
/*   double eps = 1.0e-5; */
/*   double h = eps; */
/*   vector<double> g; */
/*   for (std::size_t i = 0; i < x.size(); i++) { */
/*     if (std::abs(x[i]) > eps) */
/*       h = eps * std::abs(x[i]); */
/*     else */
/*       h = eps; */

/*     x[i] = x[i] + h; */
/*     double ff = getValue(x); */
/*     x[i] = x[i] - h; */

/*     x[i] = x[i] - h; */
/*     double fb = getValue(x); */
/*     x[i] = x[i] + h; */
/*     double df = (ff - fb)/(2.0 * h); */
/*     g.push_back(df); */
/*   } */
/*   return g; */
/* } */

vector<double> BoundaryEdgeReshapeObjFunc :: getInitialGuess()
{
  return getParamCoords();
}

void BoundaryEdgeReshapeObjFunc :: setNodes(const vector<double> &x)
{
  vector<apf::Vector3> en;
  vector<apf::Vector3> fn;
  vector<apf::Vector3> tn;
  vector<apf::Vector3> nod = convertParamCoordsToNodeVector(x);

  int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (int i = 0; i < nEN; i++)
    en.push_back(nod[i]);

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  int kk = 0;
  int nFN = mesh->getShape()->countNodesOn(mesh->TRIANGLE);

  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    for (int j = 0; j < nFN; j++)
      fn.push_back(nod[nEN+i*nFN+j]);
  }

  if (d > 2 && P > 3) {
   for (std::size_t i = nEN+nFN; i <nod.size(); i++)
     tn.push_back(nod[i]);
  }

  updateNodes(en, fn, tn, false);
}

void BoundaryEdgeReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ien, ifn, itn, true);
}


void BoundaryEdgeReshapeObjFunc :: getInitEdgeN()
{
  apf::Vector3 intEdgeX;
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));

  for (int i = 0; i < numENodes; i++) {
    mesh->getPoint(edge, i, intEdgeX);
    ien.push_back(intEdgeX);
  }
}

void BoundaryEdgeReshapeObjFunc :: getInitFaceN()
{
  apf::Adjacent adjF;
  apf::Vector3 intFaceX;
  apf::Vector3 ipFN;
  mesh->getAdjacent(edge, 2, adjF);
  int numFNodes = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    for (int j = 0; j < numFNodes; j++) {
      mesh->getPoint(adjF[i], j, intFaceX);
      ifn.push_back(intFaceX);
      //ipFN = getInterpolatingPointOnFace(mesh, adjF[i], P, j);
      //itpfn.push_back(ipFN);
    }
  }
}

void BoundaryEdgeReshapeObjFunc :: getInitTetN()
{
  apf::Adjacent adjT;
  apf::Vector3 intTetX;
  mesh->getAdjacent(edge, 3, adjT);
  int numTNodes = mesh->getShape()->countNodesOn(mesh->TET);
  for (std::size_t i = 0; i < adjT.getSize(); i++) {
    for (int j = 0; j < numTNodes; j++) {
      mesh->getPoint(adjT[i], j, intTetX);
      itn.push_back(intTetX);
    }
  }
}

vector<double> BoundaryEdgeReshapeObjFunc :: getParamCoords()
{
  apf::Vector3 xi;
  apf::Vector3 param;
  std::vector<double> xp;

  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (int i = 0; i < numENodes; i++) {
    getBezierNodeXi(mesh->getType(edge), P, i, xi);
    transferParametricOnEdgeSplit(mesh, edge, 0.5*(xi[0]+1.0), param);
    for (int j = 0; j < 3; j++) {
      xp.push_back(param[j]);
    }
  }

  apf::Vector3 xif;
  apf::Vector3 paramf;
  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
    if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      for (int j = 0; j < numFNodes; j++) {
	getBezierNodeXi(mesh->getType(adjF[i]), P, j, xif);
	transferParametricOnTriSplit(mesh, adjF[i], xif, paramf);
	for (int k = 0; k < 3; k++) {
	  xp.push_back(paramf[k]);
	}
      }
    }
    else {
      apf::Vector3 intFN;
      for (int j = 0; j < numFNodes; j++) {
	mesh->getPoint(adjF[i], j, intFN);
	for (int k = 0; k < 3; k++)
	  xp.push_back(intFN[k]);
      }
    }
  }
  return xp;
}

vector<apf::Vector3> BoundaryEdgeReshapeObjFunc :: convertParamCoordsToNodeVector(
    const vector<double> &x)
{
  apf::ModelEntity* me = mesh->toModel(edge);
  std::vector<apf::Vector3> edn;
  std::vector<apf::Vector3> vn = convertXtoNodeVector(x);
  int nENodes = mesh->getShape()->countNodesOn(mesh->EDGE);
  for (int i = 0; i < nENodes; i++) {
    apf::Vector3 coorde;
    mesh->snapToModel(me, vn[i], coorde);
    edn.push_back(coorde);
  }

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  int numFNTotal = 0;
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    int nFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));

    if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      apf::ModelEntity* mef = mesh->toModel(adjF[i]);
      for (int j = 0; j < nFNodes; j++) {
	apf::Vector3 coordf;
	mesh->snapToModel(mef, vn[nENodes+numFNTotal], coordf);
        edn.push_back(coordf);
        numFNTotal++;
      }
    }
    else {
      for (int j = 0; j < nFNodes; j++) {
	edn.push_back(vn[nENodes+numFNTotal]);
	numFNTotal++;
      }
    }
  }
  return edn;
}

vector<apf::Vector3> BoundaryEdgeReshapeObjFunc :: convertXtoNodeVector(
    const std::vector<double> &x)
{
  std::vector<apf::Vector3> a;
  apf::Vector3 v;

  if (d == 3) {
    std::size_t num = x.size()/d;
    //int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (std::size_t i = 0; i < num; i++) {
      v = {x[d*i], x[d*i + 1], x[d*i + 2]};
      a.push_back(v);
    }
  }
  return a;
}



void BoundaryEdgeReshapeObjFunc :: blendTris(
    const vector<apf::Vector3> &egn,
    vector<apf::Vector3> &faceNodes)
{
  apf::Vector3 xi;
  apf::Adjacent adjF;
  apf::Adjacent adjE;

  mesh->getAdjacent(edge, 2, adjF);

  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      mesh->getAdjacent(adjF[i], 1, adjE);
      int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));

      for (std::size_t j = 0; j < adjE.getSize(); j++) {
	if (adjE[j] == edge) {
	  int jj = 1;

	  if ( j == 0)
	    jj = 2;
	  else if ( j == 1)
	    jj = 0;
	  else
	    jj = 1;

	  for (int k = 0; k < numFNodes; k++) {
	    getBezierNodeXi(mesh->TRIANGLE, P, k, xi);
	    for (std::size_t ii = 0; ii < egn.size(); ii++) {
	      double factor =
	      	0.5 * (xi[j]/(1-xi[jj])) * binomial(P, ii+1) * intpow(1-xi[jj], ii+1) * intpow(xi[jj], P-ii-1)
	      + 0.5 * (xi[jj]/(1-xi[j])) * binomial(P, ii+1) * intpow(xi[j], ii+1) * intpow(1-xi[j], P-ii-1);
	      faceNodes[numFNodes*i+k] = faceNodes[numFNodes*i+k] + (egn[ii] - ien[ii])*factor;
	    }
	  }
	}
      }
    }
  }
}


vector<apf::Vector3> BoundaryEdgeReshapeObjFunc :: getFaceControlPointsFromInterpolatingPoints(
    apf::MeshEntity* face,
    const vector<apf::Vector3> &faceInterpolatingP)
{
  vector<apf::Vector3> faceControlP;
  apf::Vector3 xi;
  apf::NewArray<apf::Vector3> allCntrlP;
  apf::Element* Fa = apf::createElement(mesh->getCoordinateField(), face);
  apf::getVectorNodes(Fa, allCntrlP);

  int n = mesh->getShape()->countNodesOn(mesh->getType(face));
  mth::Matrix<double> A(n, n);
  mth::Matrix<double> Ainv(n, n);
  apf::NewArray<apf::Vector3> rhs(n);
  int j = 0;

  for (int i = 0; i < n; i++) {
    getBezierNodeXi(apf::Mesh::TRIANGLE, P, i, xi);
    rhs[i].zero();
    for (int ii = 0; ii < P+1; ii++) {
      for (int jj = 0; jj < P+1-ii; jj++) {
	if (ii == 0 || jj == 0 || (ii+jj == P)) {
	  double bFactor = trinomial(P, ii, jj) * Bijk(ii, jj, P-ii-jj, 1.-xi[0]-xi[1], xi[0], xi[1]);
          rhs[i] += allCntrlP[getTriNodeIndex(P, ii, jj)] * bFactor;
	}
	else {
	  j = getTriNodeIndex(P, ii, jj) - 3*P; //3P is the total number of nodes on all edges
	  A(i, j) = trinomial(P, ii, jj) * Bijk(ii, jj, P-ii-jj, 1.-xi[0]-xi[1], xi[0], xi[1]);
	}
      }
    }
    rhs[i] = faceInterpolatingP[i] - rhs[i];
  }
  apf::destroyElement(Fa);

  if (n > 1)
    invertMatrixWithPLU(n, A, Ainv);
  else
    Ainv(0,0) = 1./A(0,0);

  for (int i = 0; i < n; i++) {
    apf::Vector3 fcp(0., 0., 0.);
    for (int j = 0; j < n; j++)
      fcp += rhs[j]*Ainv(i, j);
    faceControlP.push_back(fcp);
  }
  return faceControlP;
}

void BoundaryEdgeReshapeObjFunc :: updateNodes(
    const vector<apf::Vector3> &ed,
    const vector<apf::Vector3> &fa,
    const vector<apf::Vector3> &te,
    bool isInitialX)
{
  if (!isInitialX) {
    apf::NewArray<apf::Vector3> eIntpCords;
    apf::Element* Ed = apf::createElement(mesh->getCoordinateField(), edge);
    apf::getVectorNodes(Ed, eIntpCords);
    apf::NewArray<double> trsCoff;

    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
    apf::NewArray<apf::Vector3> contP(nEN);

    for (int i = 0; i < nEN; i++) {
      for (int j = 0; j < 3; j++)
	eIntpCords[2+i][j] = ed[i][j];
    }

    getBezierTransformationCoefficients(P, mesh->getType(edge), trsCoff);
    crv::convertInterpolationPoints(nEN+2, nEN, eIntpCords, trsCoff, contP);

    for (int i = 0; i < nEN; i++)
      mesh->setPoint(edge, i, contP[i]);

    apf::destroyElement(Ed);
  }
  else {
    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (int i = 0; i < nEN; i++)
      mesh->setPoint(edge, i, ed[i]);
  }

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  std::vector<apf::Vector3> fNd;

  if (!isInitialX) {
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
      fNd.clear();

      for (int j = 0; j < numFNodes; j++)
	fNd.push_back(fa[numFNodes*i + j]);

      if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
	std::vector<apf::Vector3> fCN = getFaceControlPointsFromInterpolatingPoints(adjF[i], fNd);

	for (int j = 0; j < numFNodes; j++)
	  mesh->setPoint(adjF[i], j, fCN[j]);
      }
      else {
	for (int j = 0; j < numFNodes; j++)
	  mesh->setPoint(adjF[i], j, fa[numFNodes*i + j]);
      }
    }
  }
  else {
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
      for (int j = 0; j < numFNodes; j++)
	mesh->setPoint(adjF[i], j, fa[numFNodes*i+j]);
    }
  }

  if (d == 3 && P > 3) {
    apf::Adjacent adjT;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      int numTNodes = mesh->getShape()->countNodesOn(mesh->getType(adjT[i]));
      for (int j =0; j < numTNodes; j++)
        mesh->setPoint(adjT[i], j, te[numTNodes*i + j]);
    }
  }

  apf::synchronize(mesh->getCoordinateField());
}





// Face Reshape Objective Function
int FaceReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

double FaceReshapeObjFunc :: getValue(const std::vector<double> &x)
{
  setNodes(x);

  double sum = 0.0;

  if (d == 3) {
    apf::Adjacent adjT;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(face, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      sum = sum + computeFValNIJKL(mesh, adjT[i]);
    }
    restoreInitialNodes();
  }
  return sum;
}

/* vector<double> FaceReshapeObjFunc :: getGrad(const vector<double> &_x) */
/* { */
/*   vector<double> x; */
/*   for (int i = 0; i < _x.size(); i++) { */
/*     x[i] = _x[i]; */
/*   } */

/*   //double fold = getValue(x); */
/*   // TODO : Why the following value */
/*   double eps = 1.0e-4; */
/*   double h = eps; */
/*   std::vector<double> g; */
/*   double xmx = x[0]; */
/*   double xmn = x[0]; */
/*   double ymx = x[1]; */
/*   double ymn = x[1]; */
/*   double zmx = x[2]; */
/*   double zmn = x[2]; */
/*   double df = 0.0, dff = 0.0, dbb = 0.0; */

/*   for (std::size_t i = 0; i < x.size(); i+=3) { */
/*     if (x[i] >= xmx) xmx = x[i]; */
/*     if (x[i] <= xmn) xmn = x[i]; */
/*   } */

/*   for (std::size_t i = 1; i < x.size(); i+=3) { */
/*     if (x[i] >= ymx) ymx = x[i]; */
/*     if (x[i] <= ymn) ymn = x[i]; */
/*   } */

/*   for (std::size_t i = 2; i < x.size(); i+=3) { */
/*     if (x[i] >= zmx) zmx = x[i]; */
/*     if (x[i] <= zmn) zmn = x[i]; */
/*   } */

/*   double delx = std::abs(xmx - xmn); */
/*   double dely = std::abs(ymx - ymn); */
/*   double delz = std::abs(zmx - zmn); */
/*   double delta = 1.0; */

/*   for (std::size_t i = 0; i < x.size(); i++) { */
/*     if (i % 3 == 0) delta = delx; */
/*     if (i % 3 == 1) delta = dely; */
/*     if (i % 3 == 2) delta = delz; */

/*     if (delta < eps) */
/*       h = eps * std::abs(x[i]); */
/*     else */
/*       h = eps * delta; */

/*     x[i] = x[i] + h; */
/*     double ff = getValue(x); */
/*     x[i] = x[i] - h; */

/*     x[i] = x[i] - h; */
/*     double fb = getValue(x); */
/*     x[i] = x[i] + h; */

/*     df = (ff - fb)/(2.0 * h); */
/*     g.push_back(df); */
/*   } */
/*   return g; */
/* } */

vector<double> FaceReshapeObjFunc :: getInitialGuess()
{
  return convertNodeVectorToX(ifn);
}

void FaceReshapeObjFunc :: setNodes(const vector<double> &x)
{
  vector<apf::Vector3> fn;// (ifn.begin(), ifn.end());
  vector<apf::Vector3> tn;// (itn.begin(), itn.end());
  vector<apf::Vector3> nod = convertXtoNodeVector(x);
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
void FaceReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ifn, itn);
}

void FaceReshapeObjFunc :: getInitFaceN()
{
  apf::Vector3 intFaceX;
  int numFNodes = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  for (int j = 0; j < numFNodes; j++) {
      mesh->getPoint(face, j, intFaceX);
      ifn.push_back(intFaceX);
  }
}

void FaceReshapeObjFunc :: getInitTetN()
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

vector<double> FaceReshapeObjFunc :: convertNodeVectorToX(const vector<apf::Vector3> &fn)
{
  vector<double> x0;

  for (std::size_t i = 0; i < fn.size(); i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(fn[i][j]);

  return x0;
}

vector<apf::Vector3> FaceReshapeObjFunc :: convertXtoNodeVector(const vector<double> &x)
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

void FaceReshapeObjFunc :: updateNodes(
    const vector<apf::Vector3> &fa,
    const vector<apf::Vector3> &te)
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

}
