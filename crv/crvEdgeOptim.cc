#include "crvEdgeOptim.h"
#include "LBFGS.h"
#include "crv.h"
#include "gmi.h"
#include "apfMDS.h"
#include "apfNumbering.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvMath.h"
#include <iostream>
#include "apfMatrix.h"

/* static int global_counter = 0; */
/* static apf::MeshEntity* tetra[100]; */
/* static int number = 0; */

static void printTetNumber(apf::Mesh2* m, apf::MeshEntity* e)
{
  apf::Numbering* n = m->findNumbering("debug_num_tet");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, e, 0, 0);
  std::cout<<" TET:: "<< num <<std::endl;
}
static void printInvalidities(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* edge, int nat)
{
  apf::Numbering* n = m->findNumbering("debug_num_edge");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, edge, 0, 0);
  std::cout<<"at edge "<< num << std::endl;
  for (int i = 0; i < nat; i++) {
    std::vector<int> ai = crv::getAllInvalidities(m, e[i]);
    for (std::size_t j = 0; j < ai.size(); j++) {
      printf("%d ", ai[j]);
    }
    printf("\n");
  }
}


static void makeMultipleEntityMesh(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* edge, const char* prefix, int nat)
{
  apf::Numbering* n = m->findNumbering("debug_num_edge");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, edge, 0, 0);
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
     
 // std::stringstream ss2;
 // ss2 << "straight_sided_" << count;
 // apf::writeVtkFiles(ss2.str().c_str(), outMesh);
 
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

static void makeIndividualTetsFromFacesOrEdges(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* edge, const char* prefix, int nat)
{
  apf::Numbering* n = m->findNumbering("debug_num_edge");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, edge, 0, 0);
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

int CrvEdgeReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

void CrvEdgeReshapeObjFunc :: getInitEdgeN()
{
  apf::Vector3 intEdgeX;
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));

  for (int i = 0; i < numENodes; i++) {
    mesh->getPoint(edge, i, intEdgeX);
    ien.push_back(intEdgeX);
  }
}

void CrvEdgeReshapeObjFunc :: getInitFaceN()
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

void CrvEdgeReshapeObjFunc :: getInitTetN()
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

std::vector<double> CrvEdgeReshapeObjFunc :: getInitialGuess()
{
  return convertNodeVectorToX(ien, ifn, itn);
}

std::vector<double> CrvEdgeReshapeObjFunc :: convertNodeVectorToX(std::vector<apf::Vector3> en, std::vector<apf::Vector3> fn, std::vector<apf::Vector3> tn)
{
  std::vector<double> x0;
  for (int i = 0; i < P-1; i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(en[i][j]);

  for (std::size_t i = 0; i < fn.size(); i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(fn[i][j]);
  
  if (d > 2 && P > 3) {
    for (std::size_t i = 0; i < tn.size(); i++)
      for (int j = 0; j < 3; j++)
      	x0.push_back(tn[i][j]);
  }
  return x0;
}

std::vector<apf::Vector3> CrvEdgeReshapeObjFunc :: convertXtoNodeVector(const std::vector<double> &x)
{
  std::vector<apf::Vector3> a;
  apf::Vector3 v;
  if (d == 2) {
    std::size_t num = x.size()/d;
    //int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    //check later for 2D case:: x should not include z coordinate in optim search
    for (std::size_t i = 0; i < num; i++) {
      v = {x[d*i], x[d*i + 1], 0.0};
      a.push_back(v);
    }
  }
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

void CrvEdgeReshapeObjFunc :: blendTris(const std::vector<apf::Vector3> &egn, std::vector<apf::Vector3> &faceNodes)
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
	    double factor = 0.5 * (xi[j]/(1-xi[jj])) * binomial(P, ii+1) * intpow(1-xi[jj], ii+1) * intpow(xi[jj], P-ii-1)
	      + 0.5 * (xi[jj]/(1-xi[j])) * binomial(P, ii+1) * intpow(xi[j], ii+1) * intpow(1-xi[j], P-ii-1);
	    faceNodes[numFNodes*i+k] = faceNodes[numFNodes*i+k] + (egn[ii] - cien[ii])*factor;
	  }
	}
      }
    }
  }
}

void CrvEdgeReshapeObjFunc :: updateNodes(std::vector<apf::Vector3> ed, std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te)
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


void CrvEdgeReshapeObjFunc :: setNodes(std::vector<double> &x)
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

static double getAr(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2)
{
  p1 = p1 - p0;
  p2 = p2 - p0;
  double area = (apf::cross(p1, p2)).getLength();
  return (area/2.0);
}

std::vector<double> CrvEdgeReshapeObjFunc :: getVolume()
{
  if (d == 3) {
    apf::Adjacent adjT;
    apf::Matrix3x3 m;
    mesh->getAdjacent(edge, 3, adjT);
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
  if (d == 2) {
    apf::Adjacent adjF;
    apf::Matrix3x3 m;
    mesh->getAdjacent(edge, 2, adjF);
    apf::Vector3 point;
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      apf::Adjacent adjV;
      mesh->getAdjacent(adjF[i], 0, adjV);
      for (std::size_t j = 0; j < adjV.getSize(); j++) {
	mesh->getPoint(adjV[j], 0, point);
	for (int k = 0; k < 3; k++) {
	  if (k < 2) m[j][k] = point[k];
	  else m[j][2] = 1.0;
	}
      }
      double v = getDeterminant(m);
      vol.push_back(v);
    }
  }
  return vol;
}

double CrvEdgeReshapeObjFunc :: computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm)
{
  int weight = 1;
  double sumf = 0;
  if (d == 3) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  for (int L = 0; L <= d*(P-1); L++) {
	    if ((I == J && J == K && I == 0) || (J == K && K == L && J == 0) || (I == K && K == L && I == 0) || (I == J && J == L && I == 0))
	      weight = 14;
	    else if ((I == J && I == 0) || (I == K && I == 0) || (I == L && I == 0) || (J == K && J == 0) || (J == L && J == 0) || (K == L && K == 0))
	      weight = 4;
	    else
	      weight = 1;
	    if (I + J + K + L == d*(P-1)) {
	      double f = Nijkl(nodes,P,I,J,K)/(6.0*volm) - 1.0;
	      //std::cout<<"["<<I<<","<<J<<","<<K<<","<<L<<"]   "<<f<<std::endl;
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
	  if ((I == J && I == 0) || (J == K && J == 0) || (I == K && I == 0))
	    weight = 2;
	  else
	    weight = 1;
	  if (I + J + K == d*(P-1)) {
	    double f = Nijk(nodes,P,I,J)/(4.0*volm) - 1.0;
	    sumf = sumf + weight*f*f;
	  }
	}
      }
    }
  }
  return sumf;
}

void CrvEdgeReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ien, ifn, itn);
}

double CrvEdgeReshapeObjFunc :: getValue(std::vector<double> &x)
{
  setNodes(x);

  double sum = 0.0;
  if (d == 2) {
    apf::Adjacent adjF;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 2, adjF);
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      apf::Element* el = apf::createElement(mesh->getCoordinateField(), adjF[i]);
      apf::getVectorNodes(el, allNodes);
      sum = sum + computeFValOfElement(allNodes, vol[i]);
      apf::destroyElement(el);
    }

restoreInitialNodes();
  }

  if (d == 3) {
    apf::Adjacent adjT;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      apf::Element* el = apf::createElement(mesh->getCoordinateField(), adjT[i]);
      apf::getVectorNodes(el, allNodes);
      sum = sum + computeFValOfElement(allNodes, vol[i]);
      apf::destroyElement(el);
    }
/*
    apf::NewArray<apf::Vector3> eCP;
    apf::Element* Edl = apf::createElement(mesh->getCoordinateField(), edge);
    apf::getVectorNodes(Edl, eCP);

    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));

    double xr = 1.0;
    double xir = 1.0;
    //double alpha = 1.0;
    double beta = 0.0;
    double ad = 0.0;
    apf::Vector3 x1, x0;
    apf::Vector3 xi2, xi1, xi0;

    for (int i = 0; i < nEN; i++) {
      getBezierNodeXi(mesh->getType(edge), P, i, xi1);
      if (i > 0 && i < nEN - 1) {
        getBezierNodeXi(mesh->getType(edge), P, i+1, xi2);
        getBezierNodeXi(mesh->getType(edge), P, i-1, xi0);
        x1 = eCP[2+i+1] - eCP[2+i];
        x0 = eCP[2+i] - eCP[2+i-1];
        xir = (xi2[0] - xi1[0])/(xi1[0] - xi0[0]);
      }
      else if ( i == 0) {
        getBezierNodeXi(mesh->getType(edge), P, i+1, xi2);
        x1 = eCP[2+i+1] - eCP[2+i];
        x0 = eCP[2+i] - eCP[0];
        xir = (xi2[0] - xi1[0])/(xi1[0] + 1); // parent coordinate[-1,1]
      }
      else {
        getBezierNodeXi(mesh->getType(edge), P, i-1, xi0);
        x1 = eCP[1] - eCP[2+i];
        x0 = eCP[2+i] - eCP[2+i-1];
        xir = (1 - xi1[0])/(xi1[0] - xi0[0]);
      }

     // if (0.5*(xi1[0]+1.0) < 1.0 - 0.5*(xi1[0]+1.0))
     //   alpha = 0.5*(xi1[0]+1.0);
     // else
     //   alpha = 1.0 - 0.5*(xi1[0]+1.0);

      xr = (x1.getLength()/x0.getLength());
      ad = (xr/xir - 1);//(alpha*alpha);
      beta = beta + ad*ad;

      //sum = sum + ad*ad;
    }

    apf::destroyElement(Edl);
    */
/*
    apf::Adjacent adjF;
    mesh->getAdjacent(edge, 2, adjF);
    apf::NewArray<apf::Vector3> fCP;
    apf::Vector3 xif;
    double a[3] = {1.0, 1.0, 1.0};
    double b[3] = {1.0, 1.0, 1.0};
    double aratio = 0.0;
    //double wfactor = 1.0;
    double gamma = 0.0;

    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      apf::Element* Fal = apf::createElement(mesh->getCoordinateField(), adjF[i]);
      apf::getVectorNodes(Fal, fCP);
      int nFN = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));

      int vN[3] = {getTriNodeIndex(P, P, 0), getTriNodeIndex(P, 0, P), getTriNodeIndex(P, 0, 0)};

      double triAphys = getAr(fCP[0], fCP[1], fCP[2]);
      apf::Vector3 prt0 = {1, 0, 0};
      apf::Vector3 prt1 = {0, 1, 0};
      apf::Vector3 prt2 = {0, 0, 1};
      double triAparnt = getAr(prt0, prt1, prt2);

      for (int j = 0; j < nFN; j++) {
        getBezierNodeXi(mesh->getType(adjF[i]), P, j, xif);
        apf::Vector3 xifm = {1.0-xif[0]-xif[1], xif[0], xif[1]};
        a[0] = getAr(fCP[3*P + j], fCP[vN[0]], fCP[vN[1]]);
        b[0] = getAr(xifm, prt0, prt1);
        a[1] = getAr(fCP[3*P + j], fCP[vN[1]], fCP[vN[2]]);
        b[1] = getAr(xifm, prt1, prt2);
        a[2] = getAr(fCP[3*P + j], fCP[vN[2]], fCP[vN[0]]);
        b[2] = getAr(xifm, prt2, prt1);

        //for (int jj = 0; jj < 3; jj++) {
        //  if (xifm[i] < wfactor)
        //    wfactor = xifm[i];
        //}

        for (int k = 0; k < 3; k++) {
          aratio = (a[k]*triAparnt/(b[k]*triAphys) - 1.0); //(wfactor*wfactor);
          gamma = gamma + aratio*aratio;

          //sum = sum + aratio *aratio;
        }
      }
      //sum = sum*(1 + gamma);
      //std::cout<<"sum and ratio------ "<< sum <<"  "<<aratio-3.0<<std::endl; 
      apf::destroyElement(Fal);
    }
*/
    //sum = sum*(1 + beta);// + gamma);

    restoreInitialNodes();
  }
  return sum;
}

std::vector<double> CrvEdgeReshapeObjFunc :: getGrad(std::vector<double> &x)
{
  //double fold = getValue(x);
  double eps = 1.0e-5;
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
/*
    dff = (ff - fold)/h;
    dbb = (fold - fb)/h;
    df = (dff + dbb)/(2.0);
*/
    df = (ff - fb)/(2.0 * h);
    g.push_back(df);
  }
  return g;
}

void CrvEdgeOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvEdgeOptim :: setTol(double tolerance)
{
  tol = tolerance;
}

bool CrvEdgeOptim :: run(int &invaliditySize)
{
  std::vector<int> sizeHolder;
  apf::MeshEntity* adj_array[99];
  apf::Adjacent adj;
  mesh->getAdjacent(edge, 3, adj);
  for (int i = 0; i < adj.getSize(); i++) {
    adj_array[i] = adj[i];
    std::vector<int> ai = crv::getAllInvalidities(mesh, adj[i]);
    sizeHolder.push_back(ai.size());
  }

  std::vector<int> ai = crv::getAllInvalidities(mesh, tet);
  //makeMultipleEntityMesh(mesh, adj_array, edge, "before_cavity_of_edge_", adj.getSize());
  //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "before_cavity_indv_tet_of_edge_", adj.getSize());
  printTetNumber(mesh, tet);
  printInvalidities(mesh, adj_array, edge, adj.getSize());
  CrvEdgeReshapeObjFunc *objF = new CrvEdgeReshapeObjFunc(mesh, edge, tet);
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  LBFGS *l = new LBFGS(tol, iter, x0, objF);

  apf::MeshEntity* ed[6];

  for (std::size_t i = 0; i < adj.getSize(); i++) {
    mesh->getDownward(adj[i], 1, ed);
    int edgeIndex = apf::findIn(ed, 6, edge);
    printf("reshape tried on %d edge; ", edgeIndex);
    printTetNumber(mesh, adj[i]);
  }

  bool hasDecreased = false;

  if (l->run()) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);


    for (int i = 0; i < adj.getSize(); i++) {
      std::vector<int> aiNew = crv::getAllInvalidities(mesh, adj[i]);
      invaliditySize = invaliditySize + aiNew.size();
      hasDecreased = hasDecreased || (aiNew.size() > sizeHolder[i]);
    }

    if (hasDecreased == false) {
      //makeMultipleEntityMesh(mesh, adj_array, edge, "after_cavity_of_edge_", adj.getSize());
      //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "after_cavity_indv_tet_of_edge_", adj.getSize());
      printInvalidities(mesh, adj_array, edge, adj.getSize());
      std::cout<<"--------------------------------------"<<std::endl;
      return true;
    }
    else {
      objF->restoreInitialNodes();
      printInvalidities(mesh, adj_array, edge, adj.getSize());
      std::cout<<"Size DID NOT decrease"<<std::endl;
      std::cout<<"--------------------------------------"<<std::endl;
      return false;
    }
  }
  else {
    //finalX = l->currentX; 
    //objF->setNodes(finalX); 
    //makeMultipleEntityMesh(mesh, adj_array, edge, "after_cavity_of_edge_", adj.getSize());
    std::cout<<"*****Edge Optim FAILURE" <<std::endl;
    std::cout<<"--------------------------------------"<<std::endl;
    return false;
  }
}

}
