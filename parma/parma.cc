#include <PCU.h>
#include "parma.h"
#include <parma_dcpart.h>
#include <limits>

namespace {
  int numSharedFaces(apf::Mesh* m) {
    apf::MeshIterator *it = m->begin(m->getDimension()-1);
    apf::MeshEntity* e;
    int cnt = 0;
    while( (e = m->iterate(it)) )
      if( m->isShared(e) )
        cnt++;
    m->end(it);
    return cnt;
  }
  int numOwnedVtx(apf::Mesh* m) {
    apf::MeshIterator *it = m->begin(0);
    apf::MeshEntity* e;
    int cnt = 0;
    while( (e = m->iterate(it)) )
      if( m->isShared(e) && m->isOwned(e) )
        cnt++;
    m->end(it);
    return cnt;
  }
}

void Parma_GetEntImbalance(apf::Mesh* mesh, double (*entImb)[4]) {
   double tot[4];
   for(int i=0; i<= mesh->getDimension(); i++)
      tot[i] = (*entImb)[i] = mesh->count(i);
   if( mesh->getDimension() != 3 )
      tot[3] = (*entImb)[3] = 0;
   PCU_Add_Doubles(tot, 4);
   PCU_Max_Doubles(*entImb, 4);
   for(int i=0; i<4; i++)
      (*entImb)[i] /= (tot[i]/PCU_Comm_Peers());
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
  loc = neighbors.size();
  max = loc;
  PCU_Max_Ints(&max,1);
  double total = loc;
  PCU_Add_Doubles(&total,1);
  avg = total / PCU_Comm_Peers();
}

long Parma_GetNumBdryVtx(apf::Mesh* m) {
  long cnt = numOwnedVtx(m);
  PCU_Add_Longs(&cnt, 1);
  return cnt;
}

void Parma_GetDisconnectedStats(apf::Mesh* m, int& max, double& avg, int& loc) {
  dcPart dc(m);
  int tot = max = loc = dc.numDisconnectedComps();
  PCU_Debug_Print("getDisStats dc parts %d\n", loc);
  PCU_Max_Ints(&max, 1);
  PCU_Add_Ints(&tot, 1);
  avg = static_cast<double>(tot)/PCU_Comm_Peers();
}

void Parma_ProcessDisconnectedParts(apf::Mesh* m) {
  dcPart dc(m);
  dc.fix();
}

void Parma_PrintPtnStats(apf::Mesh* m, std::string key) {
  int maxDc = 0;
  double avgDc = 0;
  int locDc = 0;
  Parma_GetDisconnectedStats(m, maxDc, avgDc, locDc);

  int maxNb = 0;
  double avgNb = 0;
  int locNb = 0;
  Parma_GetNeighborStats(m, maxNb, avgNb, locNb);

  const long bdryVtx = Parma_GetNumBdryVtx(m);

  int surf = numSharedFaces(m);
  int vol = m->count(m->getDimension());
  double minSurfToVol, maxSurfToVol, avgSurfToVol;
  minSurfToVol =  maxSurfToVol =  avgSurfToVol = surf/(double)vol;
  PCU_Min_Doubles(&minSurfToVol, 1);
  PCU_Max_Doubles(&maxSurfToVol, 1);
  PCU_Add_Doubles(&avgSurfToVol, 1);
  avgSurfToVol /= PCU_Comm_Peers();

  int empty = (m->count(m->getDimension()) == 0 ) ? 1 : 0;
  PCU_Add_Ints(&empty, 1);

  double imb[4] = {0, 0, 0, 0};
  Parma_GetEntImbalance(m, &imb);

  if( 0 == PCU_Comm_Self() ) {
    fprintf(stdout, "STATUS %s disconnected <max avg> %d %.3f\n",
        key.c_str(), maxDc, avgDc);
    fprintf(stdout, "STATUS %s neighbors <max avg> %d %.3f\n",
        key.c_str(), maxNb, avgNb);
    fprintf(stdout, "STATUS %s empty parts %d\n",
        key.c_str(), empty);
    fprintf(stdout, "STATUS %s number of shared vtx %ld\n",
        key.c_str(), bdryVtx);
    fprintf(stdout, "STATUS %s sharedFacesToElements <max min avg> "
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

