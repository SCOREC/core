#include "parma.h"
#include <parma_diffmc.h>
#include <parma_dcpart.h>
#include <parma_imbinfo.h>
#include <limits>
#include <PCU.h>

int Parma_Run(apf::Mesh* mesh, apf::MeshTag* weight, const double maxImb) {
  double entImb[4];
  Parma_GetWeightedEntImbalance(mesh, weight, &entImb);
  int priority[4] = {0, 0, 0, 1};
  int ierr = Parma_RunWeightedPtnImprovement(mesh, weight, &priority);
  return ierr;
}

int Parma_RunPtnImprovement(apf::Mesh* mesh, int (*priority)[4], 
                const double maxImb, const int dbgLvl, const int maxItr){
   Parma parma(mesh);
   int ierr = parma.run(priority, dbgLvl, maxItr, maxImb);
   return ierr;
}

int Parma_RunWeightedPtnImprovement(apf::Mesh* mesh, 
                apf::MeshTag* weight, int (*priority)[4], 
                const double maxImb, const int dbgLvl, 
                const int maxItr) {
   Parma parma(mesh, weight);
   int ierr = parma.run(priority, dbgLvl, maxItr, maxImb);
   return ierr;
}

void Parma_GetEntImbalance(apf::Mesh* mesh, double (*entImb)[4]) {
   imbInfo imb;
   apf::MeshTag* weight = mesh->createDoubleTag("weight",1);
   imb.get(mesh,weight);
   for(int i=0; i<= mesh->getDimension(); i++) 
      (*entImb)[i] = imb.maxImb[i];
   mesh->destroyTag(weight); 
}

void Parma_GetWeightedEntImbalance(apf::Mesh* mesh, apf::MeshTag* weight, double (*entImb)[4]) {
   imbInfo imb;
   imb.get(mesh,weight);
   for(int i=0; i<4; i++) 
      (*entImb)[i] = imb.maxImb[i];
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
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity* e;
  long cnt = 0;
  while ((e = m->iterate(it))) 
    if (m->isShared(e) && m->isOwned(e))
      cnt++;
  m->end(it);
  PCU_Add_Longs(&cnt, 1);
  return cnt;
}

void Parma_GetDisconnectedStats(apf::Mesh* m, int& max, double& avg, int& loc) {
  dcPart dc(m);
  loc = dc.numDisconnectedComps();
  int tot = max = loc;
  PCU_Max_Ints(&max, 1);
  PCU_Add_Ints(&tot, 1);
  avg = static_cast<double>(tot)/PCU_Comm_Peers();
}

void Parma_ProcessDisconnectedParts(apf::Mesh* m) {
  dcPart dc(m);
  int numTotDc = dc.numDisconnectedComps();
  PCU_Add_Ints(&numTotDc, 1);
  if ( numTotDc > 0 ) {
    if( 0 == PCU_Comm_Self() ) 
      fprintf(stderr, "PARMA_STATUS initial number of disconnected components %d\n", numTotDc);
    dc.fix();
    numTotDc = dc.numDisconnectedComps();
    PCU_Add_Ints(&numTotDc, 1);
    if( 0 == PCU_Comm_Self() ) 
      fprintf(stderr, "PARMA_STATUS after fix() number of disconnected components %d\n", numTotDc);
  }
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

  int empty = (m->count(m->getDimension()) == 0 ) ? 1 : 0;
  PCU_Add_Ints(&empty, 1);

  double imb[4] = {0, 0, 0, 0};
  Parma_GetEntImbalance(m, &imb);

  if( 0 == PCU_Comm_Self() ) {
    fprintf(stdout, "STATUS %s disconnected <max avg> %d %.3f\n", key.c_str(), maxDc, avgDc);
    fprintf(stdout, "STATUS %s neighbors <max avg> %d %.3f\n", key.c_str(), maxNb, avgNb);
    fprintf(stdout, "STATUS %s empty parts %d\n", key.c_str(), empty);
    fprintf(stdout, "STATUS %s number of shared vtx %ld\n", key.c_str(), bdryVtx);
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

