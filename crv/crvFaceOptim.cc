#include "crvFaceOptim.h"
#include "LBFGS.h"
#include "crv.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvMath.h"
#include <iostream>
#include "apfMatrix.h"
#include "apfNumbering.h"
#include "gmi.h"
#include "apfMDS.h"

static void printTetNumber(apf::Mesh2* m, apf::MeshEntity* e)
{
 // printf("\n");
 // return;
  apf::Numbering* n = m->findNumbering("debug_num_tet");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, e, 0, 0);
  std::cout<<" TET:: "<< num <<std::endl;
}

static void printInvalidities(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* face, int nat)
{

//  return;
  apf::Numbering* n = m->findNumbering("debug_num_face");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, face, 0, 0);
  int tag = m->getModelTag(m->toModel(face));
  std::cout<<"at face "<< num <<" tag: "<<tag<< std::endl;

  for (int i = 0; i < nat; i++) {
    std::vector<int> ai = crv::getAllInvalidities(m, e[i]);
    for (std::size_t j = 0; j < ai.size(); j++) {
      printf("%d ", ai[j]);
    }
    printf("\n");
  }
}

static void makeMultipleEntityMesh(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* face, const char* prefix, int nat)
{
  return;
  apf::Numbering* n = m->findNumbering("debug_num_face");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, face, 0, 0);
  int dim = 0;
  if (m->getType(e[0]) == apf::Mesh::TRIANGLE)
    dim = 2;
  else if (m->getType(e[0]) == apf::Mesh::TET)
    dim = 3;
  else
    PCU_ALWAYS_ASSERT(0);

  gmi_model* g = gmi_load(".null");
  apf::Mesh2* outMesh = apf::makeEmptyMdsMesh(g, dim, false);

  apf::MeshEntity* newEnt[99];


  for (int ii = 0; ii < nat; ii++) {
    // Verts
    apf::MeshEntity* vs[12];
    apf::MeshEntity* newVs[12];
    int nv = m->getDownward(e[ii], 0, vs);
    for(int i = 0; i < nv; ++i)
    {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[i], 0, p);
      newVs[i] = outMesh->createVertex(0, p, param);
    }

    // Edges
    apf::MeshEntity* es[12];
    apf::MeshEntity* newEs[12];
    int ne = m->getDownward(e[ii], 1, es);
    for(int i = 0; i < ne; ++i)
    {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* new_evs[2];
      m->getDownward(es[i], 0, evs);
      for (int j = 0; j < 2; j++) {
        new_evs[j] = newVs[apf::findIn(vs, nv, evs[j])];
      }

      newEs[i] = outMesh->createEntity(apf::Mesh::EDGE, 0, new_evs);
    }

    // Faces
    apf::MeshEntity* fs[12];
    apf::MeshEntity* newFs[12];
    int nf = m->getDownward(e[ii], 2, fs);
    for(int i = 0; i < nf; ++i)
    {
      apf::MeshEntity* fes[3];
      apf::MeshEntity* new_fes[3];
      m->getDownward(fs[i], 1, fes);
      for (int j = 0; j < 3; j++) {
        new_fes[j] = newEs[apf::findIn(es, ne, fes[j])];
      }
      newFs[i] = outMesh->createEntity(apf::Mesh::TRIANGLE, 0, new_fes);
    }

    // Regions
    apf::MeshEntity* tet;
    if (dim == 3) {
      tet = outMesh->createEntity(apf::Mesh::TET, 0, newFs);
    }

    if (dim == 2)
      newEnt[ii] = newFs[0];
    else
      newEnt[ii] = tet;

    PCU_ALWAYS_ASSERT(m->getType(e[ii]) == outMesh->getType(newEnt[ii]));
    //printf("HERE 02\n")liver tetrahedral element
    outMesh->acceptChanges();
  }
  apf::changeMeshShape(outMesh, crv::getBezier(3), true);
  outMesh->acceptChanges();

  for (int ii = 0; ii < nat; ii++) {
    for (int d = 1; d <= dim; d++)
    {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* eds[12];
      int counter = m->getDownward(e[ii], d, eds);
      apf::MeshEntity* new_eds[12];
      outMesh->getDownward(newEnt[ii], d, new_eds);
      int non = outMesh->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for(int n = 0; n < counter; ++n) {
        for(int i = 0; i < non; ++i) {
        apf::Vector3 p;
        m->getPoint(eds[n], i, p);
        outMesh->setPoint(new_eds[n], i, p);
        }
      }
    }
    outMesh->acceptChanges();
  }

  std::stringstream ss;
  ss << prefix<< num;
  crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TET, 8, ss.str().c_str());
  crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TRIANGLE, 8, ss.str().c_str());
  crv::writeCurvedWireFrame(outMesh, 8, ss.str().c_str());

  outMesh->destroyNative();
  apf::destroyMesh(outMesh);
}

static void visualizeAllFacesOfTet(apf::Mesh2* m, apf::MeshEntity* e, int count, const char* prefix)
{
  int dim = 0;
  if (m->getType(e) == apf::Mesh::TRIANGLE) {
    std::cout<<"the entity is not a TET"<<std::endl;
  }
  else if (m->getType(e) == apf::Mesh::TET)
    dim = 3;
  else
    PCU_ALWAYS_ASSERT(0);

  gmi_model* g = gmi_load(".null");

  apf::Mesh2* outMesh[4];
  for (int i = 0; i < 4; i++) {
    outMesh[i] = apf::makeEmptyMdsMesh(g, 2, false);
  }

  apf::MeshEntity* face[4];
  apf::MeshEntity* newface[4];
  int nf = m->getDownward(e, 2, face);
 for (int i = 0; i < nf; i++) {

    //Verts
    apf::MeshEntity* vs[3];
    apf::MeshEntity* newVs[3];
    int nv = m->getDownward(face[i], 0, vs);
    for (int j = 0; j < nv; j++) {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[j], 0, p);
      newVs[j] = outMesh[i]->createVertex(0, p, param);
    }
    outMesh[i]->acceptChanges();

    //Edges
    apf::MeshEntity* es[3];
    apf::MeshEntity* newEs[3];
    int ne = m->getDownward(face[i], 1, es);
    for (int j = 0; j < ne; j++) {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* newEvs[2];
      m->getDownward(es[j], 0, evs);
      for (int k = 0; k < 2; k++) {
        int kk = apf::findIn(vs,nv, evs[k]);
        newEvs[k] = newVs[kk];
      }

      newEs[j] = outMesh[i]->createEntity(apf::Mesh::EDGE, 0, newEvs);
    }

    //Faces
    apf::MeshEntity* fes[3];
    apf::MeshEntity* newFes[3];
    m->getDownward(face[i], 1, fes);
    for (int j = 0; j < 3; j++)
      newFes[j] = newEs[apf::findIn(es, ne, fes[j])];

    newface[i] = outMesh[i]->createEntity(apf::Mesh::TRIANGLE, 0, newFes);

    PCU_ALWAYS_ASSERT(m->getType(face[i]) == outMesh[i]->getType(newface[i]));
    outMesh[i]->acceptChanges();

    apf::changeMeshShape(outMesh[i], crv::getBezier(3), true);
    outMesh[i]->acceptChanges();

    for (int d = 1; d < dim; d++) {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* ent[10];
      int counter = m->getDownward(face[i], d, ent);
      apf::MeshEntity* newent[10];
      outMesh[i]->getDownward(newface[i], d, newent);
      int non = outMesh[i]->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for (int n = 0; n < counter; n++) {
        for (int j = 0; j < non; j++) {
          apf::Vector3 p;
          m->getPoint(ent[n], j, p);
          outMesh[i]->setPoint(newent[n], j, p);
        }
      }
    }

    outMesh[i]->acceptChanges();

    std::stringstream ss;
    ss << prefix << "_Face_"<< i;
    crv::writeCurvedVtuFiles(outMesh[i], apf::Mesh::TRIANGLE, 40, ss.str().c_str());
    crv::writeCurvedWireFrame(outMesh[i], 50, ss.str().c_str());
  }
}

static void makeIndividualTetsFromFacesOrEdges(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* face, const char* prefix, int nat)
{
  apf::Numbering* n = m->findNumbering("debug_num_face");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, face, 0, 0);
  int dim = 0;
  if (m->getType(e[0]) == apf::Mesh::TRIANGLE)
    dim = 2;
  else if (m->getType(e[0]) == apf::Mesh::TET)
    dim = 3;
  else
    PCU_ALWAYS_ASSERT(0);

  gmi_model* g = gmi_load(".null");

  apf::Mesh2* outMesh[nat];
  for (int i = 0; i < nat; i++) {
    outMesh[i] = apf::makeEmptyMdsMesh(g, dim, false);
  }

  apf::MeshEntity* newEnt[99];

  for (int ii = 0; ii < nat; ii++) {
    // Verts
    apf::MeshEntity* vs[12];
    apf::MeshEntity* newVs[12];
    int nv = m->getDownward(e[ii], 0, vs);
    for(int i = 0; i < nv; ++i)
    {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[i], 0, p);
      newVs[i] = outMesh[ii]->createVertex(0, p, param);
    }

    // Edges
    apf::MeshEntity* es[12];
    apf::MeshEntity* newEs[12];
    int ne = m->getDownward(e[ii], 1, es);
    for(int i = 0; i < ne; ++i)
    {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* new_evs[2];
      m->getDownward(es[i], 0, evs);
      for (int j = 0; j < 2; j++) {
        new_evs[j] = newVs[apf::findIn(vs, nv, evs[j])];
      }

      newEs[i] = outMesh[ii]->createEntity(apf::Mesh::EDGE, 0, new_evs);
    }
   // Faces
    apf::MeshEntity* fs[12];
    apf::MeshEntity* newFs[12];
    int nf = m->getDownward(e[ii], 2, fs);
    for(int i = 0; i < nf; ++i)
    {
      apf::MeshEntity* fes[3];
      apf::MeshEntity* new_fes[3];
      m->getDownward(fs[i], 1, fes);
      for (int j = 0; j < 3; j++) {
        new_fes[j] = newEs[apf::findIn(es, ne, fes[j])];
      }
      newFs[i] = outMesh[ii]->createEntity(apf::Mesh::TRIANGLE, 0, new_fes);
    }

    // Regions
    apf::MeshEntity* tet;
    if (dim == 3) {
      tet = outMesh[ii]->createEntity(apf::Mesh::TET, 0, newFs);
    }

    if (dim == 2)
      newEnt[ii] = newFs[0];
    else
      newEnt[ii] = tet;

    PCU_ALWAYS_ASSERT(m->getType(e[ii]) == outMesh[ii]->getType(newEnt[ii]));
    outMesh[ii]->acceptChanges();

    apf::changeMeshShape(outMesh[ii], crv::getBezier(3), true);
    outMesh[ii]->acceptChanges();


    for (int d = 1; d <= dim; d++)
    {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* eds[12];
      int counter = m->getDownward(e[ii], d, eds);
      apf::MeshEntity* new_eds[12];
      outMesh[ii]->getDownward(newEnt[ii], d, new_eds);
      int non = outMesh[ii]->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for(int n = 0; n < counter; ++n) {
        for(int i = 0; i < non; ++i) {
        apf::Vector3 p;
        m->getPoint(eds[n], i, p);
        outMesh[ii]->setPoint(new_eds[n], i, p);
        }
      }
    }
    outMesh[ii]->acceptChanges();

    std::stringstream ss;
    ss << prefix<< num << "_TET_"<<ii;
    crv::writeCurvedVtuFiles(outMesh[ii], apf::Mesh::TET, 8, ss.str().c_str());
    crv::writeCurvedVtuFiles(outMesh[ii], apf::Mesh::TRIANGLE, 8, ss.str().c_str());
    crv::writeCurvedWireFrame(outMesh[ii], 8, ss.str().c_str());

    //outMesh[ii]->destroyNative();
    //apf::destroyMesh(outMesh[ii]);

    visualizeAllFacesOfTet(m, e[ii], ii, ss.str().c_str());
  }
}




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
	      weight = 6;
	    else if ((I == J && I == 0) || (I == K && I == 0) || (I == L && I == 0) || (J == K && J == 0) || (J == L && J == 0) || (K == L && K == 0))
	      weight = 3;
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

void CrvFaceOptim :: setTol(double tolerance)
{
  tol = tolerance;
}

bool CrvFaceOptim :: run(int &invaliditySize)
{
  std::vector<int> sizeHolder;
  apf::MeshEntity* adj_array[99];
  apf::Adjacent adj;
  mesh->getAdjacent(face, 3, adj);
  int thisTetSize = 0;
  for (int i = 0; i < adj.getSize(); i++) {
    adj_array[i] = adj[i];
    //std::vector<int> ai = crv::getAllInvalidities(mesh, adj[i]);
    //if (adj[i] == tet) thisTetSize = ai.size();
    //sizeHolder.push_back(ai.size());   
  }

  std::vector<int> ai = crv::getAllInvalidities(mesh, tet);
  invaliditySize = ai.size();
  //makeMultipleEntityMesh(mesh, adj_array, face, "before_cavity_of_face_", adj.getSize());
  //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, face, "before_cavity_indv_tet_of_face_", adj.getSize());
  printTetNumber(mesh, tet);
  printInvalidities(mesh, adj_array, face, adj.getSize());
  CrvFaceReshapeObjFunc *objF = new CrvFaceReshapeObjFunc(mesh, face, tet);
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  LBFGS *l = new LBFGS(tol, iter, x0, objF);

  apf::MeshEntity* fc[4];
  int thisTETnum = 0;

  for (std::size_t i = 0; i < adj.getSize(); i++) {
    if (adj[i] == tet) thisTETnum = 1;
    else thisTETnum = 0;
    mesh->getDownward(adj[i], 2, fc);
    int faceIndex = apf::findIn(fc, 4, face);
    printf("reshape tried on %d face, TET %d ",faceIndex, thisTETnum);
    printTetNumber(mesh, adj[i]);
  }

  //bool hasDecreased = false;
  //invaliditySize = 0;

  //if (l->run() && thisTetSize > 0) {
  if (l->run() && invaliditySize > 0) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);

    std::vector<int> aiNew = crv::getAllInvalidities(mesh, tet);
    /*
    for (int i = 0; i < adj.getSize(); i++) {
      std::vector<int> aiNew = crv::getAllInvalidities(mesh, adj[i]);
      invaliditySize = invaliditySize + aiNew.size();
      hasDecreased = hasDecreased || (aiNew.size() > sizeHolder[i]);
    }
*/
   // if (hasDecreased == false ) {
    if (aiNew.size() <= invaliditySize) {
      invaliditySize = aiNew.size();
      //makeMultipleEntityMesh(mesh, adj_array, face, "after_cavity_of_face_", adj.getSize());
      //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, face, "after_cavity_indv_tet_of_face_", adj.getSize());
      printInvalidities(mesh, adj_array, face, adj.getSize());
      std::cout<<"----------------------------------------------------"<<std::endl;

      return true;
    }
    else {
      objF->restoreInitialNodes();
      printInvalidities(mesh, adj_array, face, adj.getSize());
      std::cout<<"Size didNOT decrease"<<std::endl;
      std::cout<<"----------------------------------------------------"<<std::endl;
      return false;
    }
  }
  else {
    if (thisTetSize == 0) {
      std::cout<<" No Optimization tried"<<std::endl;
      std::cout<<"-------------------------------------------"<<std::endl;
      return false;
    }
    std::cout<<"*****FaceOptim FAILURE"<<std::endl;
    return false;
  }
}

}
