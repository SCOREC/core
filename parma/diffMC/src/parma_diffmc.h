#ifndef PARMA_DIFFMC_H_
#define PARMA_DIFFMC_H_

#include "apf.h"
#include "apfMesh.h"
#include "parma_diffmcversion.h"
#include "parma_priority.h"
#include "parma_partinfo.h"
#include <list>

#define ParMA_Embed_Version() \
        {static const char *junk="ParMA Diffusive Multi Criteria Improvement version " \
        PARMADIFFMC_MAJOR_VERSION "." PARMADIFFMC_MINOR_VERSION "." PARMADIFFMC_PATCH_VERSION ; } ;

#define ParMA_Print_Version() \
        {fprintf(stdout, "ParMA Diffusive Multi Criteria Improvement version %s\n", \
        PARMADIFFMC_MAJOR_VERSION "." PARMADIFFMC_MINOR_VERSION "." PARMADIFFMC_PATCH_VERSION); } ;

typedef apf::DynamicArray<apf::MeshEntity*> eArr;
typedef std::list<apf::MeshEntity*> eList;

class Parma {
   public:  
      Parma(apf::Mesh*& m);
      Parma(apf::Mesh*& m, apf::MeshTag* tag);
      ~Parma();
      int run(int (*priority)[4], const int dbgLvl, 
          const int maxIter, const double maxImb);
   protected: 
      apf::Mesh* mesh;
   private:
      Parma();
      void init(apf::Mesh*& m);
      void improveBalance(const priorityList &pl, const int plIdx);
      double improveEntBalance(const priorityList &pl, const int plIdx, 
            apf::Array<double,4> avgW, const int itr);

      double tagElmsForMigr(partInfo& part, const double maxW, apf::Migration* plan);
      double tagElmsForMigr(partInfo& part, const double maxW, 
          const int maxFace, apf::Migration* plan);

      double tagElmsForMigr_Vtx(partInfo& part, const double maxW, apf::Migration* plan);
      double tagSmallCavitiesForMigr(partInfo& part, const double maxW, 
          const int maxAdjElm, apf::Migration* plan);
      double tagPtnMdlEdgeCavities(partInfo& part, const double maxW, 
          const int maxAdjElm, apf::Migration* plan);

      bool inputsValid(int (*priority)[4], int dbgLvl, int maxIter, 
          double maxImb);

      apf::MeshTag* wtag;
      apf::MeshTag* vtag;
      bool userWeights;

      int maxIter;
      int dbgLvl;
      double maxImb;
};

#endif
