/****************************************************************************** 

  (c) 2013-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "SimParasolidKrnl.h"
#include "SimAcisKrnl.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "SimMessages.h"
#include "ModelTypes.h"
#include "MeshTypes.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

using std::string;

void messageHandler(int type, const char *msg);
void progressHandler(const char *what, int level, int startVal, 
                     int endVal, int currentVal, void *);

/* forward declare */
pGModel GM_translateAcisModel(pGModel, double tolerance);
pGModel GM_translateParasolidModel(pGModel, double tolerance);

int main(int argc, char* argv[])
{
  string acisExt = ".sat";
  string paraExt = ".xmt_txt";

  if( argc != 2 ) {
     fprintf(stderr, "Usage: %s <acis %s or parasolid %s model file>\n", 
         argv[0], acisExt.c_str(), paraExt.c_str());
     return 0;
  }

  pGModel model;

  try {  
    Sim_readLicenseFile(0);
    Sim_logOn("relations.log");
    SimModel_start();
    SimParasolid_start(1);
    SimAcis_start(1);

    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setCallback(progress, progressHandler);
    
    // Load the model as an assembly model
    string mdlName = argv[1];
    bool isAcis = (mdlName.find(acisExt) != std::string::npos);
    bool isParasolid = (mdlName.find(paraExt) != std::string::npos);
    if( isAcis ) {
      pAcisNativeModel acisNModel;
      acisNModel = AcisNM_createFromFile(mdlName.c_str(),0);
      model = GM_createFromNativeModel(acisNModel,progress);
      NM_release(acisNModel);
    } else if ( isParasolid ) {
      pParasolidNativeModel paraNModel;
      paraNModel = ParasolidNM_createFromFile(mdlName.c_str(),0);
      model = GM_createFromNativeModel(paraNModel,progress);
      NM_release(paraNModel);
    } else {
      fprintf(stderr, "ERROR model extension is neither %s or %s... exiting\n", 
          acisExt.c_str(), paraExt.c_str());
      exit(EXIT_FAILURE);
    }
   
    // check that the input model is topologically valid
    pPList modelErrors = PList_new();
    if(!GM_isValid(model,0,modelErrors)) {
      printf("Input model not valid\n");
      printf("Number of errors returned: %d\n",PList_size(modelErrors));
      GM_release(model);
      return 1;
    }
    PList_delete(modelErrors);

    // translate the model
    pGModel simmodel;
    if( isAcis )
      simmodel = GM_translateAcisModel(model, 1e-04);
    else if( isParasolid )
      simmodel = GM_translateParasolidModel(model, 1e-04);

    GM_write(model, "model.smd", 0, progress);
    GM_write(simmodel, "translated-model.smd", 0, progress);
   
    GM_release(simmodel);
    GM_release(model);
    Progress_delete(progress);
    
    SimParasolid_stop(1);
    SimAcis_stop(1);
    SimModel_stop();
    Sim_logOff();
    Sim_unregisterAllKeys();
    
  } catch (pSimError err) {
    printf("Simmetrix error caught:\n");
    printf("  Error code: %d\n",SimError_code(err));
    printf("  Error string: %s\n",SimError_toString(err));
    SimError_delete(err);
    return 1;
  } catch (...) {
    printf("Unhandled exception caught\n");
    return 1;
  }
  return 0;
}

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    printf("Info: %s\n",msg);
    break;
  case Sim_DebugMsg:
    printf("Debug: %s\n",msg);
    break;
  case Sim_WarningMsg:
    printf("Warning: %s\n",msg);
    break;
  case Sim_ErrorMsg:
    printf("Error: %s\n",msg);
    break;
  }
  return;
}

void progressHandler(const char *what, int level, int startVal, 
                     int endVal, int currentVal, void *)
{
  printf("Progress: %s, level: %d, startVal: %d, endVal: %d, currentVal: %d\n",
         what,level,startVal,endVal,currentVal);
  return;
}
