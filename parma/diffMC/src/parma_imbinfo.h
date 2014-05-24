#ifndef PARMA_IMBINFO_H_
#define PARMA_IMBINFO_H_

#include "apf.h"
#include "apfMesh.h"

class imbInfo {
   public:
      double maxImb[4];
      double maxW[4];
      double totW[4];
      apf::Array<double, 4> minAvgW;
      apf::Array<double, 4> curAvgW;
      double w[4];
   
      imbInfo();
      void init();
      void get(apf::Mesh* m, apf::MeshTag* wtag);
      void print();
   private:
      void getWeight(apf::Mesh* m, apf::MeshTag* wtag);
      void getWeightedImbalance();
      bool printedHeader;
};

#endif

