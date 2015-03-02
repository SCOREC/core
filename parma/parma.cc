#include <PCU.h>
#include "parma.h"
#include "diffMC/maximalIndependentSet/mis.h"
#include <parma_dcpart.h>
#include <limits>
#include <assert.h>

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
  double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
    assert(m->hasTag(e,w));
    double weight;
    m->getDoubleTag(e,w,&weight);
    return weight;
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

  long totEnt[4] = {0,0,0,0};
  int minEnt[4] = {0,0,0,0}, maxEnt[4] = {0,0,0,0};
  double avgEnt[4] = {0,0,0,0};
  for( int d=0; d<=m->getDimension(); d++)
    entStats(m, d, totEnt[d], minEnt[d], maxEnt[d], avgEnt[d]);

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
  Parma_GetEntImbalance(m, &imb);

  if (fine) {
    fprintf(stdout, "FINE STATUS %s <Partid vtx rgn dc nb "
                    "owned_bdry shared_bdry model_bdry shSidesToElm > "
                    " %d %lu %lu %d %d %d %d %d %.3f\n",
      key.c_str(), PCU_Comm_Self()+1, m->count(0), m->count(m->getDimension()),
      locDc, locNb, locV[0], locV[1], locV[2], surf/TO_DBL(vol));
    PCU_Barrier();
  }

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

    const char* orders[4] = {"vtx","edge","face","rgn"};
    for( int d=0; d<=m->getDimension(); d++)
      fprintf(stdout, "STATUS %s %s <tot max min avg> "
          "%ld %d %d %.3f\n",
          key.c_str(), orders[d], totEnt[d], maxEnt[d], minEnt[d], avgEnt[d]);
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
