#include "PCU.h"
#include "parma_imbinfo.h"
#include "parma_commons.h"
#include <assert.h>
#include <vector>
#include <limits>
#include "mpi.h"

using apf::MeshIterator;

using parmaCommons::isLess;
using parmaCommons::status;

imbInfo::imbInfo() {
  init();
}

void imbInfo::init() {
  for(int i=0; i<4; i++) {
    maxImb[i]=0;
    maxW[i]=0;
    totW[i]=0;
    minAvgW[i]=std::numeric_limits<double>::max();
    curAvgW[i]=0;
    w[i]=0;
  }
  printedHeader = false;
}

void imbInfo::get(apf::Mesh* m, apf::MeshTag* wtag) {
  getWeight(m, wtag);
  getWeightedImbalance();
}

void imbInfo::print() {
  if ( 0 == PCU_Comm_Self() ) {
    if( ! printedHeader ) {
      status("imbInfo maxVtxImb maxEdgeImb maxFaceImb maxRgnImb "
          "minVtxAvg minEdgeAvg minFaceAvg minRgnAvg "
          "vtxAvg edgeAvg faceAvg rgnAvg "
          "totVtxW totEdgeW totFaceW totRgnW\n");
      printedHeader = true;
    }
    status("imbInfo %.3lf %.3lf %.3lf %.3lf "
        "%.3lf %.3lf %.3lf %.3lf "
        "%.3lf %.3lf %.3lf %.3lf "
        "%.3lf %.3lf %.3lf %.3lf\n", 
        maxImb[0], maxImb[1], maxImb[2], maxImb[3], 
        minAvgW[0], minAvgW[1], minAvgW[2], minAvgW[3], 
        curAvgW[0], curAvgW[1], curAvgW[2], curAvgW[3], 
        totW[0], totW[1], totW[2], totW[3] );
  }
}

void imbInfo::getWeight(apf::Mesh* m, apf::MeshTag* wtag) {
  apf::MeshEntity* e;
  for(int d=0; d<=m->getDimension(); d++) { 
    w[d] = 0;
    MeshIterator* itr = m->begin(d);
    while( (e = m->iterate(itr)) ) {
      double weight = 1;
      if ( m->hasTag(e, wtag) ) {
        m->getDoubleTag(e, wtag, &weight);
      }
      w[d] += weight;
    }
    m->end(itr);
  }
}

void imbInfo::getWeightedImbalance() {
  int commSz;
  PCU_Comm_Size(&commSz);

  for(int i=0; i<4; i++) 
    totW[i] = maxW[i] = w[i];
  PCU_Add_Doubles(totW, 4); 
  PCU_Max_Doubles(maxW, 4); 

  for ( int i=0; i<4; i++ ) {
    curAvgW[i] = totW[i] / (double)commSz;
    if( isLess(curAvgW[i], minAvgW[i]) ) 
      minAvgW[i] = curAvgW[i]; 
    maxImb[i] = (double) maxW[i] / (double) minAvgW[i];
  }
}



