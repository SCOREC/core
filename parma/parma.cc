#include <pcu_util.h>
#include "parma.h"
#include "diffMC/maximalIndependentSet/mis.h"
#include "diffMC/parma_commons.h"
#include "diffMC/parma_convert.h"
#include <parma_dcpart.h>
#include <lionPrint.h>
#include <limits>
#include <sstream>
#include <string>

namespace {
  typedef std::map<int,int> mii;

  int numSharedSides(apf::Mesh* m) {
    apf::MeshIterator *it = m->begin(m->getDimension()-1);
    apf::MeshEntity* e;
    int cnt = 0;
    while( (e = m->iterate(it)) )
      if( m->isShared(e) )
        cnt++;
    m->end(it);
    return cnt;
  }

  int numBdryVtx(apf::Mesh* m, bool onlyShared=false) {
    apf::MeshIterator *it = m->begin(0);
    apf::MeshEntity* e;
    int cnt = 0;
    while( (e = m->iterate(it)) )
      if( m->isShared(e) && (onlyShared || m->isOwned(e)) )
          cnt++;
    m->end(it);
    return cnt;
  }

  int numMdlBdryVtx(apf::Mesh* m) {
    const int dim = m->getDimension();
    apf::MeshIterator *it = m->begin(0);
    apf::MeshEntity* e;
    int cnt = 0;
    while( (e = m->iterate(it)) )
      if( m->getModelType(m->toModel(e)) < dim )
          cnt++;
    m->end(it);
    return cnt;
  }

  void hasEntWeight(apf::Mesh* m, apf::MeshTag* w, int (*hasWeight)[4]) {
    for(size_t i=0; i < 4; i++)
      (*hasWeight)[i] = 1;
    int dims = m->getDimension() + 1;
    for(int i=0; i < dims; i++) {
      apf::MeshIterator* it = m->begin(i);
      apf::MeshEntity* e;
      while ((e = m->iterate(it)))
        if(! m->hasTag(e,w) ) {
          (*hasWeight)[i] = 0;
          break;
        }
      m->end(it);
    }
  }

  double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
    PCU_ALWAYS_ASSERT(m->hasTag(e,w));
    double weight;
    m->getDoubleTag(e,w,&weight);
    return weight;
  }

  void getPartWeights(apf::Mesh* m, apf::MeshTag* w, double (*weight)[4]) {
    int hasWeight[4];
    hasEntWeight(m,w,&hasWeight);
    int dims = m->getDimension() + 1;
    for(int i=0; i < dims; i++) {
      (*weight)[i] = 0;
      if(hasWeight[i]) {
        apf::MeshIterator* it = m->begin(i);
        apf::MeshEntity* e;
        while ((e = m->iterate(it)))
          (*weight)[i] += getEntWeight(m, e, w);
        m->end(it);
      } else {
        (*weight)[i] += TO_DOUBLE(m->count(i));
      }
    }
  }

  void getWeightedStats(
      pcu::PCU *PCUObj, int dim, double (*loc)[4], double (*tot)[4],
      double (*min)[4], double (*max)[4], double (*avg)[4]) {
    for(int d=0; d<=dim; d++)
      (*min)[d] = (*max)[d] = (*tot)[d] = (*loc)[d];
    PCUObj->Min<double>(*min, dim+1);
    PCUObj->Max<double>(*max, dim+1);
    PCUObj->Add<double>(*tot, dim+1);
    for(int d=0; d<=dim; d++) {
      (*avg)[d] = (*tot)[d];
      (*avg)[d] /= TO_DOUBLE(PCUObj->Peers());
    }
  }

  void getStats(pcu::PCU *PCUObj, int& loc, long& tot, int& min, int& max, double& avg) {
    min = PCUObj->Min<int>(loc);
    max = PCUObj->Max<int>(loc);
    tot = PCUObj->Add<long>(TO_LONG(loc));
    avg = TO_DOUBLE(tot) / PCUObj->Peers();
  }

  using parmaCommons::status;

  void writeFineStats(apf::Mesh* m, std::string key,
      int locDc, int locNb, int* locV, int surf, double vol) {
    const char* entNames[4] = {"vtx", "edge", "face", "rgn"};
    const int dim = m->getDimension();
    std::stringstream ss;
    ss << "FINE STATUS " << key << "<Partid ";
    for(int d=0; d<=dim; d++)
      ss << entNames[d] << ' ';
    ss << "dc nb owned_bdry shared_bdry model_bdry shSidesToElm > "
       << m->getPCU()->Self()+1  << ' ';
    for(int d=0; d<=dim; d++)
      ss << m->count(d) << ' ';
    ss << m->count(m->getDimension())  << ' '
       << locDc  << ' '
       << locNb  << ' '
       << locV[0]  << ' ' <<  locV[1]  << ' ' <<  locV[2]  << ' '
       << surf/TO_DOUBLE(vol);
    std::string s = ss.str();
    lion_eprint(1, "%s\n", s.c_str());
    m->getPCU()->Barrier();
  }

  void writeWeightedEntStats(apf::Mesh* m, apf::MeshTag* w, std::string key) {
    double weight[4];
    getPartWeights(m, w, &weight);
    double minEnt[4] = {0,0,0,0}, maxEnt[4] = {0,0,0,0};
    double totEnt[4] = {0,0,0,0}, avgEnt[4] = {0,0,0,0};
    getWeightedStats(m->getPCU(), m->getDimension(), &weight, &totEnt, &minEnt, &maxEnt, &avgEnt);
    const char* orders[4] = {"vtx","edge","face","rgn"};
    if(!m->getPCU()->Self()) {
      for( int d=0; d<=m->getDimension(); d++)
        status("%s weighted %s <tot max min avg> "
            "%.1f %.1f %.1f %.3f\n",
            key.c_str(), orders[d],
            totEnt[d], maxEnt[d], minEnt[d], avgEnt[d]);
    }
  }

  void getNeighborCounts(apf::Mesh* m, mii& nborToShared) {
    apf::MeshIterator *it = m->begin(0);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::Parts sharers;
      m->getResidence(e,sharers);
      APF_ITERATE(apf::Parts, sharers, nbor)
        nborToShared[*nbor]++;
    }
    m->end(it);
  }
}

void Parma_GetEntImbalance(apf::Mesh* mesh, double (*entImb)[4]) {
  size_t dims;
  double tot[4];
  dims = TO_SIZET(mesh->getDimension()) + 1;
  for(size_t i=0; i < dims; i++)
    tot[i] = (*entImb)[i] = mesh->count(TO_INT(i));
  mesh->getPCU()->Add<double>(tot, dims);
  mesh->getPCU()->Max<double>(*entImb, dims);
  for(size_t i=0; i < dims; i++)
    (*entImb)[i] /= (tot[i]/mesh->getPCU()->Peers());
  for(size_t i=dims; i < 4; i++)
    (*entImb)[i] = 1.0;
}

void Parma_GetWeightedEntImbalance(apf::Mesh* mesh, apf::MeshTag* w,
    double (*entImb)[4]) {
  size_t dims = TO_SIZET(mesh->getDimension()) + 1;
  getPartWeights(mesh, w, entImb);
  double tot[4] = {0,0,0,0};
  for(size_t i=0; i < dims; i++)
    tot[i] = (*entImb)[i];
  mesh->getPCU()->Add<double>(tot, TO_SIZET(dims));
  mesh->getPCU()->Max<double>(*entImb, TO_SIZET(dims));
  for(size_t i=0; i < dims; i++)
    (*entImb)[i] /= (tot[i]/mesh->getPCU()->Peers());
  for(size_t i=dims; i < 4; i++)
    (*entImb)[i] = 1.0;
}

double Parma_GetWeightedEntImbalance(apf::Mesh* m, apf::MeshTag* w,
    int dim) {
    PCU_ALWAYS_ASSERT(dim >= 0 && dim <= 3);
    apf::MeshIterator* it = m->begin(dim);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      sum += getEntWeight(m, e, w);
    m->end(it);
   double tot = m->getPCU()->Add<double>(sum);
   double max = m->getPCU()->Max<double>(sum);
   return max/(tot/m->getPCU()->Peers());
}

void Parma_GetNeighborStats(apf::Mesh* m, int& max, int& numMaxParts,
    double& avg, int& loc) {
  mii nborToShared;
  getNeighborCounts(m,nborToShared);
  loc = TO_INT(nborToShared.size())-1;
  max = m->getPCU()->Max<int>(loc);
  avg = TO_DOUBLE(m->getPCU()->Add<int>(loc)) / m->getPCU()->Peers();
  numMaxParts = m->getPCU()->Add<int>( (loc==max) );
}

void Parma_WriteSmallNeighbors(apf::Mesh* m, int small, const char* prefix) {
  mii nborToShared;
  getNeighborCounts(m,nborToShared);
  int* smallCnt = new int[small];
  for(int i=0; i<small; i++) smallCnt[i] = 0;
  APF_ITERATE(mii, nborToShared, nbor)
    for(int i=0; i<small; i++)
      if( nbor->second == i+1 )
        smallCnt[i]++;
  m->getPCU()->Add<int>(smallCnt,small);
  if( !m->getPCU()->Self() ) {
    std::stringstream ss;
    for(int i=0; i<small; i++)
      ss << i+1 << ":" << smallCnt[i] << " ";
    std::string s = ss.str();
    status("%s small neighbor counts %s\n", prefix, s.c_str());
  }
  delete [] smallCnt;
}

int Parma_GetSmallestSideMaxNeighborParts(apf::Mesh* m) {
  mii nborToShared;
  getNeighborCounts(m,nborToShared);
  int loc = TO_INT(nborToShared.size())-1;
  int max = m->getPCU()->Max<int>(loc);
  int smallest = INT_MAX;
  if( loc == max ) {
    APF_ITERATE(mii, nborToShared, nbor)
      if( nbor->second < smallest )
        smallest = nbor->second;
  }
  return m->getPCU()->Min<int>(smallest);
}

void Parma_GetOwnedBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg) {
  loc = numBdryVtx(m);
  getStats(m->getPCU(), loc, tot, min, max, avg);
}

void Parma_GetSharedBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg) {
  bool onlyShared = true;
  loc = numBdryVtx(m,onlyShared);
  getStats(m->getPCU(), loc, tot, min, max, avg);
}

void Parma_GetMdlBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg) {
  loc = numMdlBdryVtx(m);
  getStats(m->getPCU(), loc, tot, min, max, avg);
}

void Parma_GetDisconnectedStats(apf::Mesh* m, int& max, double& avg, int& loc) {
  dcPart dc(m);
  loc = TO_INT(dc.getNumDcComps());
  max = m->getPCU()->Max<int>(loc);
  avg = TO_DOUBLE( m->getPCU()->Add<int>(loc) ) / m->getPCU()->Peers();
}

void Parma_ProcessDisconnectedParts(apf::Mesh* m) {
  dcPartFixer dcf(m);
}

void Parma_PrintPtnStats(apf::Mesh* m, std::string key, bool fine) {
  apf::MeshTag* w = m->createDoubleTag("parma_ent_weights", 1);
  int dims = m->getDimension() + 1;
  double entWeight=1;
  for(int i=0; i < dims; i++) {
    apf::MeshIterator* it = m->begin(i);
    apf::MeshEntity* e;
    while ((e = m->iterate(it)))
      m->setDoubleTag(e,w,&entWeight);
    m->end(it);
  }
  Parma_PrintWeightedPtnStats(m,w,key,fine);
  for(int i=0; i < dims; i++)
    apf::removeTagFromDimension(m,w,i);
  m->destroyTag(w);
}

void Parma_PrintWeightedPtnStats(apf::Mesh* m, apf::MeshTag* w, std::string key, bool fine) {
  m->getPCU()->DebugPrint("%s vtx %lu\n", key.c_str(), m->count(0));
  m->getPCU()->DebugPrint("%s edge %lu\n", key.c_str(), m->count(1));
  m->getPCU()->DebugPrint("%s face %lu\n", key.c_str(), m->count(2));
  if( m->getDimension() == 3 )
    m->getPCU()->DebugPrint("%s rgn %lu\n", key.c_str(), m->count(3));

  int maxDc = 0;
  double avgDc = 0;
  int locDc = 0;
  Parma_GetDisconnectedStats(m, maxDc, avgDc, locDc);
  m->getPCU()->DebugPrint("%s dc %d\n", key.c_str(), locDc);

  int maxNb = 0, maxNbParts = 0;
  double avgNb = 0;
  int locNb = 0;
  Parma_GetNeighborStats(m, maxNb, maxNbParts, avgNb, locNb);
  int smallSideMaxNbPart = Parma_GetSmallestSideMaxNeighborParts(m);
  m->getPCU()->DebugPrint("%s neighbors %d\n", key.c_str(), locNb);

  int locV[3], minV[3], maxV[3];
  long totV[3];
  double avgV[3];
  Parma_GetOwnedBdryVtxStats(m, locV[0], totV[0], minV[0], maxV[0], avgV[0]);
  Parma_GetSharedBdryVtxStats(m, locV[1], totV[1], minV[1], maxV[1], avgV[1]);
  Parma_GetMdlBdryVtxStats(m, locV[2], totV[2], minV[2], maxV[2], avgV[2]);
  m->getPCU()->DebugPrint("%s ownedBdryVtx %d\n", key.c_str(), locV[0]);
  m->getPCU()->DebugPrint("%s sharedBdryVtx %d\n", key.c_str(), locV[1]);
  m->getPCU()->DebugPrint("%s mdlBdryVtx %d\n", key.c_str(), locV[2]);

  int surf = numSharedSides(m);
  double vol = TO_DOUBLE( m->count(m->getDimension()) );
  double surfToVol = surf/vol;
  double minSurfToVol = m->getPCU()->Min<double>(surfToVol);
  double maxSurfToVol = m->getPCU()->Max<double>(surfToVol);
  double avgSurfToVol = m->getPCU()->Add<double>(surfToVol) / m->getPCU()->Peers();
  m->getPCU()->DebugPrint("%s sharedSidesToElements %.3f\n", key.c_str(), surfToVol);

  int empty = (m->count(m->getDimension()) == 0 ) ? 1 : 0;
  empty = m->getPCU()->Add<int>(empty);

  double imb[4] = {0, 0, 0, 0};
  Parma_GetWeightedEntImbalance(m,w,&imb);

  if (fine)
    writeFineStats(m, key, locDc, locNb, locV, surf, vol);

  m->getPCU()->DebugPrint("%s vtxAdjacentNeighbors ", key.c_str());
  apf::Parts peers;
  apf::getPeers(m,0,peers);
  APF_ITERATE(apf::Parts,peers,p)
    m->getPCU()->DebugPrint("%d ", *p);
  m->getPCU()->DebugPrint("\n");

  if( 0 == m->getPCU()->Self() ) {
    status("%s disconnected <max avg> %d %.3f\n",
        key.c_str(), maxDc, avgDc);
    status("%s neighbors <max avg> %d %.3f\n",
        key.c_str(), maxNb, avgNb);
    status("%s smallest side of max neighbor part %d\n",
        key.c_str(), smallSideMaxNbPart);
    status("%s num parts with max neighbors %d\n",
        key.c_str(), maxNbParts);
    status("%s empty parts %d\n",
        key.c_str(), empty);
  }
  const int smallestSide = 10;
  Parma_WriteSmallNeighbors(m, smallestSide, key.c_str());
  writeWeightedEntStats(m,w,key);
  if( 0 == m->getPCU()->Self() ) {
    status("%s owned bdry vtx <tot max min avg> "
        "%ld %d %d %.3f\n",
        key.c_str(), totV[0], maxV[0], minV[0], avgV[0]);
    status("%s shared bdry vtx <tot max min avg> "
        "%ld %d %d %.3f\n",
        key.c_str(), totV[1], maxV[1], minV[1], avgV[1]);
    status("%s model bdry vtx <tot max min avg> "
        "%ld %d %d %.3f\n",
        key.c_str(), totV[2], maxV[2], minV[2], avgV[2]);
    status("%s sharedSidesToElements <max min avg> "
        "%.3f %.3f %.3f\n",
        key.c_str(), maxSurfToVol, minSurfToVol, avgSurfToVol);
    status("%s entity imbalance <v e f r>: "
        "%.2f %.2f %.2f %.2f\n", key.c_str(), imb[0], imb[1], imb[2], imb[3]);
  }
}

apf::MeshTag* Parma_WeighByMemory(apf::Mesh* m) {
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  apf::MeshTag* tag = m->createDoubleTag("parma_bytes", 1);
  while ((e = m->iterate(it))) {
    double bytes = m->getElementBytes(m->getType(e));
    m->setDoubleTag(e, tag, &bytes);
  }
  m->end(it);
  return tag;
}

int Parma_MisNumbering(apf::Mesh* m, int d) {
  apf::Parts neighbors;
  apf::getPeers(m,d,neighbors);

  misLuby::partInfo part;
  part.id = m->getId();
  part.net.push_back(m->getId());
  APF_ITERATE(apf::Parts, neighbors, nItr) {
    part.adjPartIds.push_back(*nItr);
    part.net.push_back(*nItr);
  }

  unsigned seed = TO_UINT(part.id+1);
  mis_init(seed);
  int misNumber=-1;
  int iter=0;
  int misSize=0;
  while( misSize != m->getPCU()->Peers() ) {
    if( mis(part, m->getPCU(), false, true) || 1 == part.net.size() ) {
      misNumber = iter;
      part.net.clear();
      part.adjPartIds.clear();
    }
    iter++;
    misSize = (misNumber != -1);
    misSize = m->getPCU()->Add<int>(misSize);
  }
  return misNumber;
}
