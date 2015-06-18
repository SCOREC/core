#include <PCU.h>
#include "parma.h"
#include "diffMC/maximalIndependentSet/mis.h"
#include <parma_dcpart.h>
#include <limits>
#include <assert.h>
#include <sstream>
#include <string>

#define TO_SIZET(a) static_cast<size_t>(a)
#define TO_INT(a) static_cast<int>(a)
#define TO_DBL(a) static_cast<double>(a)

namespace {
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
    size_t dims = TO_SIZET(m->getDimension()) + 1;
    for(size_t i=0; i < dims; i++) {
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
    assert(m->hasTag(e,w));
    double weight;
    m->getDoubleTag(e,w,&weight);
    return weight;
  }
  void getPartWeights(apf::Mesh* m, apf::MeshTag* w, double (*weight)[4]) {
    int hasWeight[4];
    hasEntWeight(m,w,&hasWeight);
    size_t dims = TO_SIZET(m->getDimension()) + 1;
    for(size_t i=0; i < dims; i++) {
      (*weight)[i] = 0;
      if(hasWeight[i]) {
        apf::MeshIterator* it = m->begin(i);
        apf::MeshEntity* e;
        while ((e = m->iterate(it)))
          (*weight)[i] += getEntWeight(m, e, w);
        m->end(it);
      } else {
          (*weight)[i] += TO_DBL(m->count(i));
      }
    }
  }
  void getWeightedStats(double (*loc)[4], double (*tot)[4],
      double (*min)[4], double (*max)[4], double (*avg)[4]) {
    for(int d=0; d<4; d++)
      (*min)[d] = (*max)[d] = (*tot)[d] = (*loc)[d];
    PCU_Min_Doubles(*min, 4);
    PCU_Max_Doubles(*max, 4);
    PCU_Add_Doubles(*tot, 4);
    for(int d=0; d<4; d++) {
      (*avg)[d] = static_cast<double>((*tot)[d]);
      (*avg)[d] /= TO_DBL(PCU_Comm_Peers());
    }
  }
  void getStats(int& loc, long& tot, int& min, int& max, double& avg) {
    min = max = loc;
    tot = static_cast<long>(loc);
    PCU_Min_Ints(&min, 1);
    PCU_Max_Ints(&max, 1);
    PCU_Add_Longs(&tot, 1);
    avg = static_cast<double>(tot);
    avg /= TO_DBL(PCU_Comm_Peers());
  }
  void entStats(apf::Mesh* m, int dim,
      long& tot, int& min, int& max, double& avg) {
    int loc = TO_INT(m->count(dim));
    getStats(loc, tot, min, max, avg);
  }
  void writeFineStats(apf::Mesh* m, std::string key,
      int locDc, int locNb, int* locV, int surf, double vol) {
    const char* entNames[4] = {"vtx", "edge", "face", "rgn"};
    const int dim = m->getDimension();
    std::stringstream ss;
    ss << "FINE STATUS " << key << "<Partid ";
    for(int d=0; d<=dim; d++)
      ss << entNames[d] << ' ';
    ss << "dc nb owned_bdry shared_bdry model_bdry shSidesToElm > "
       << PCU_Comm_Self()+1  << ' ';
    for(int d=0; d<=dim; d++)
      ss << m->count(d) << ' ';
    ss << m->count(m->getDimension())  << ' '
       << locDc  << ' '
       << locNb  << ' '
       << locV[0]  << ' ' <<  locV[1]  << ' ' <<  locV[2]  << ' '
       << surf/TO_DBL(vol);
    std::string s = ss.str();
    fprintf(stderr, "%s\n", s.c_str());
    PCU_Barrier();
  }
}

void Parma_GetEntImbalance(apf::Mesh* mesh, double (*entImb)[4]) {
   size_t dims;
   double tot[4];
   dims = TO_SIZET(mesh->getDimension()) + 1;
   for(size_t i=0; i < dims; i++)
      tot[i] = (*entImb)[i] = mesh->count(TO_INT(i));
   PCU_Add_Doubles(tot, dims);
   PCU_Max_Doubles(*entImb, dims);
   for(size_t i=0; i < dims; i++)
      (*entImb)[i] /= (tot[i]/PCU_Comm_Peers());
   for(size_t i=dims; i < 4; i++)
      (*entImb)[i] = 1.0;
}

void Parma_PrintWeightedEntStats(apf::Mesh* m, apf::MeshTag* w, std::string key) {
  double weight[4];
  getPartWeights(m, w, &weight);
  double minEnt[4] = {0,0,0,0}, maxEnt[4] = {0,0,0,0};
  double totEnt[4] = {0,0,0,0}, avgEnt[4] = {0,0,0,0};
  getWeightedStats(&weight, &totEnt, &minEnt, &maxEnt, &avgEnt);
  const char* orders[4] = {"vtx","edge","face","rgn"};
  if(!PCU_Comm_Self()) {
    for( int d=0; d<=m->getDimension(); d++)
      fprintf(stdout, "STATUS %s weighted %s <tot max min avg> "
          "%.1f %.1f %.1f %.3f\n",
          key.c_str(), orders[d],
          totEnt[d], maxEnt[d], minEnt[d], avgEnt[d]);
  }
}

void Parma_GetWeightedEntImbalance(apf::Mesh* mesh, apf::MeshTag* w,
    double (*entImb)[4]) {
   double tot[4] = {0,0,0,0};
   size_t dims = TO_SIZET(mesh->getDimension()) + 1;
   for(size_t i=0; i < dims; i++) {
     apf::MeshIterator* it = mesh->begin(i);
     apf::MeshEntity* e;
     while ((e = mesh->iterate(it)))
       tot[i] += getEntWeight(mesh, e, w);
     (*entImb)[i] = tot[i];
     mesh->end(it);
   }
   PCU_Add_Doubles(tot, dims);
   PCU_Max_Doubles(*entImb, dims);
   for(size_t i=0; i < dims; i++)
      (*entImb)[i] /= (tot[i]/PCU_Comm_Peers());
   for(size_t i=dims; i < 4; i++)
      (*entImb)[i] = 1.0;
}

double Parma_GetWeightedEntImbalance(apf::Mesh* m, apf::MeshTag* w,
    int dim) {
    assert(dim >= 0 && dim <= 3);
    apf::MeshIterator* it = m->begin(dim);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      sum += getEntWeight(m, e, w);
    m->end(it);
   double max, tot;
   max = tot = sum;
   PCU_Add_Doubles(&tot, 1);
   PCU_Max_Doubles(&max, 1);
   return max/(tot/PCU_Comm_Peers());
}


void Parma_GetNeighborStats(apf::Mesh* m, int& max, double& avg, int& loc) {
  apf::MeshIterator *it = m->begin(0);
  apf::Parts neighbors;
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::Parts sharers;
    m->getResidence(e,sharers);
    neighbors.insert(sharers.begin(),sharers.end());
  }
  m->end(it);
  loc = static_cast<int>(neighbors.size())-1;
  max = loc;
  PCU_Max_Ints(&max,1);
  double total = loc;
  PCU_Add_Doubles(&total,1);
  avg = total / PCU_Comm_Peers();
}

void Parma_GetOwnedBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg) {
  loc = numBdryVtx(m);
  getStats(loc, tot, min, max, avg);
}

void Parma_GetSharedBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg) {
  bool onlyShared = true;
  loc = numBdryVtx(m,onlyShared);
  getStats(loc, tot, min, max, avg);
}

void Parma_GetMdlBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg) {
  loc = numMdlBdryVtx(m);
  getStats(loc, tot, min, max, avg);
}

void Parma_GetEntStats(apf::Mesh* m, int dim, long& tot, int& min, int& max,
    double& avg, int& loc) {
  assert( dim>=0 && dim<=m->getDimension() );
  entStats(m, dim, tot, min, max, avg);
  loc = TO_INT(m->count(dim));
}

void Parma_GetDisconnectedStats(apf::Mesh* m, int& max, double& avg, int& loc) {
  dcPart dc(m);
  int tot = max = loc = TO_INT(dc.getNumComps())-1;
  PCU_Max_Ints(&max, 1);
  PCU_Add_Ints(&tot, 1);
  avg = static_cast<double>(tot)/PCU_Comm_Peers();
}

void Parma_ProcessDisconnectedParts(apf::Mesh* m) {
  dcPartFixer dcf(m);
}

void Parma_PrintPtnStats(apf::Mesh* m, std::string key, bool fine) {
  apf::MeshTag* w = m->createDoubleTag("parma_ent_weights", 1);
  size_t dims = TO_SIZET(m->getDimension()) + 1;
  double entWeight=1;
  for(size_t i=0; i < dims; i++) {
    apf::MeshIterator* it = m->begin(i);
    apf::MeshEntity* e;
    while ((e = m->iterate(it)))
      m->setDoubleTag(e,w,&entWeight);
    m->end(it);
  }
  Parma_PrintWeightedPtnStats(m,w,key,fine);
  for(size_t i=0; i < dims; i++)
    apf::removeTagFromDimension(m,w,i);
  m->destroyTag(w);
}

void Parma_PrintWeightedPtnStats(apf::Mesh* m, apf::MeshTag* w, std::string key, bool fine) {
  PCU_Debug_Print("%s vtx %lu\n", key.c_str(), m->count(0));
  PCU_Debug_Print("%s edge %lu\n", key.c_str(), m->count(1));
  PCU_Debug_Print("%s face %lu\n", key.c_str(), m->count(2));
  if( m->getDimension() == 3 )
    PCU_Debug_Print("%s rgn %lu\n", key.c_str(), m->count(3));

  int maxDc = 0;
  double avgDc = 0;
  int locDc = 0;
  Parma_GetDisconnectedStats(m, maxDc, avgDc, locDc);
  PCU_Debug_Print("%s dc %d\n", key.c_str(), locDc);

  int maxNb = 0;
  double avgNb = 0;
  int locNb = 0;
  Parma_GetNeighborStats(m, maxNb, avgNb, locNb);
  PCU_Debug_Print("%s neighbors %d\n", key.c_str(), locNb);

  int locV[3], minV[3], maxV[3];
  long totV[3];
  double avgV[3];
  Parma_GetOwnedBdryVtxStats(m, locV[0], totV[0], minV[0], maxV[0], avgV[0]);
  Parma_GetSharedBdryVtxStats(m, locV[1], totV[1], minV[1], maxV[1], avgV[1]);
  Parma_GetMdlBdryVtxStats(m, locV[2], totV[2], minV[2], maxV[2], avgV[2]);
  PCU_Debug_Print("%s ownedBdryVtx %d\n", key.c_str(), locV[0]);
  PCU_Debug_Print("%s sharedBdryVtx %d\n", key.c_str(), locV[1]);
  PCU_Debug_Print("%s mdlBdryVtx %d\n", key.c_str(), locV[2]);

  int surf = numSharedSides(m);
  double vol = static_cast<double>( m->count(m->getDimension()) );
  double minSurfToVol, maxSurfToVol, avgSurfToVol, surfToVol;
  minSurfToVol =  maxSurfToVol =  avgSurfToVol = surfToVol = surf/vol;
  PCU_Min_Doubles(&minSurfToVol, 1);
  PCU_Max_Doubles(&maxSurfToVol, 1);
  PCU_Add_Doubles(&avgSurfToVol, 1);
  avgSurfToVol /= PCU_Comm_Peers();
  PCU_Debug_Print("%s sharedSidesToElements %.3f\n", key.c_str(), surfToVol);

  int empty = (m->count(m->getDimension()) == 0 ) ? 1 : 0;
  PCU_Add_Ints(&empty, 1);

  double imb[4] = {0, 0, 0, 0};
  Parma_GetWeightedEntImbalance(m,w,&imb);

  if (fine)
    writeFineStats(m, key, locDc, locNb, locV, surf, vol);

  PCU_Debug_Print("%s vtxAdjacentNeighbors ", key.c_str());
  apf::Parts peers;
  apf::getPeers(m,0,peers);
  APF_ITERATE(apf::Parts,peers,p)
    PCU_Debug_Print("%d ", *p);
  PCU_Debug_Print("\n");

  if( 0 == PCU_Comm_Self() ) {
    fprintf(stdout, "STATUS %s disconnected <max avg> %d %.3f\n",
        key.c_str(), maxDc, avgDc);
    fprintf(stdout, "STATUS %s neighbors <max avg> %d %.3f\n",
        key.c_str(), maxNb, avgNb);
    fprintf(stdout, "STATUS %s empty parts %d\n",
        key.c_str(), empty);
  }
  Parma_PrintWeightedEntStats(m,w,key);
  if( 0 == PCU_Comm_Self() ) {
    fprintf(stdout, "STATUS %s owned bdry vtx <tot max min avg> "
        "%ld %d %d %.3f\n",
        key.c_str(), totV[0], maxV[0], minV[0], avgV[0]);
    fprintf(stdout, "STATUS %s shared bdry vtx <tot max min avg> "
        "%ld %d %d %.3f\n",
        key.c_str(), totV[1], maxV[1], minV[1], avgV[1]);
    fprintf(stdout, "STATUS %s model bdry vtx <tot max min avg> "
        "%ld %d %d %.3f\n",
        key.c_str(), totV[2], maxV[2], minV[2], avgV[2]);
    fprintf(stdout, "STATUS %s sharedSidesToElements <max min avg> "
        "%.3f %.3f %.3f\n",
        key.c_str(), maxSurfToVol, minSurfToVol, avgSurfToVol);
    fprintf(stdout, "STATUS %s entity imbalance <v e f r>: "
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

  unsigned int seed = static_cast<unsigned int>(part.id+1);
  mis_init(seed);
  int misNumber=-1;
  int iter=0;
  int misSize=0;
  while( misSize != PCU_Comm_Peers() ) {
    if( mis(part, false, true) || 1 == part.net.size() ) {
      misNumber = iter;
      part.net.clear();
      part.adjPartIds.clear();
    }
    iter++;
    misSize = (misNumber != -1);
    PCU_Add_Ints(&misSize,1);
  }
  return misNumber;
}
