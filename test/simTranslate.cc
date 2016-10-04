/****************************************************************************** 

  (c) 2013-2014 Scientific Computation Research Center, 
  Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

 *******************************************************************************/
#include "SimModel.h"
#include "SimAdvModel.h"
#include "SimUtil.h"
#include "SimMessages.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "MeshSim.h"
#include "SimAttribute.h"
#include "AttributeTypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

/* cheap hackish way to get SIM_PARASOLID and SIM_ACIS */
#include "gmi_sim_config.h"
#include <gmi_sim.h>

#ifdef SIM_PARASOLID
#include "SimParasolidKrnl.h"
#endif

#ifdef SIM_ACIS
#include "SimAcisKrnl.h"
#endif

void messageHandler(int type, const char *msg);
void progressHandler(const char *what, int level, int startVal, 
    int endVal, int currentVal, void *);

void SModel_setAttManager(pModel model, pAManager attMan);
pAManager SModel_attManager(pModel model);

/* Constants (ugh -- FIXME) */
static std::string acisExt = ".sat";
static std::string paraExt = ".xmt_txt";
static std::string paraExtshort = ".x_t";

void printSimError(pSimError err)
{
  printf("Simmetrix error caught:\n");
  printf("  Error code: %d\n",SimError_code(err));
  printf("  Error string: %s\n",SimError_toString(err));
  SimError_delete(err);
}

void translateModel(std::string mdlName, pGModel* simmodel, pProgress& progress)
{
  pGModel model;
  bool isAcis = (mdlName.find(acisExt) != std::string::npos);
  bool isParasolid = (mdlName.find(paraExt) != std::string::npos);
  isParasolid |= (mdlName.find(paraExtshort) != std::string::npos);
  if( isAcis ) {
#ifdef SIM_ACIS
    pAcisNativeModel acisNModel;
    acisNModel = AcisNM_createFromFile(mdlName.c_str(),0);
    model = GM_createFromNativeModel(acisNModel,progress);
    NM_release(acisNModel);
#else
    fprintf(stderr, "model is Acis but not compiled with SIM_ACIS=ON\n");
    exit(EXIT_FAILURE);
#endif
  } else if ( isParasolid ) {
#ifdef SIM_PARASOLID
    pParasolidNativeModel paraNModel;
    paraNModel = ParasolidNM_createFromFile(mdlName.c_str(),0);
    model = GM_createFromNativeModel(paraNModel,progress);
    NM_release(paraNModel);
#else
    fprintf(stderr, "model is Parasolid but not compiled with SIM_PARASOLID=ON\n");
    exit(EXIT_FAILURE);
#endif
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
    abort();
  }
  PList_delete(modelErrors);

  // translate the model
  *simmodel = GM_translateModel(model, NULL);
  GM_release(model);
}

int main(int argc, char* argv[])
{
  pProgress progress;
  pGModel simmodel;
  pAManager attmngr;

  if( argc < 3 ) {
    printf("Usage: %s <acis %s or parasolid %s> <output .smd>\n", 
        argv[0], acisExt.c_str(), paraExt.c_str());
    printf("       to translate a native model into a GeomSim model\n");
    printf("   or: %s <acis %s or parasolid %s> <attribute .smd> <output .smd>\n", 
        argv[0], acisExt.c_str(), paraExt.c_str());
    printf("       to do the above and combine it with attributes into one file\n");
    return 0;
  }

  //Setup
  try {
    MS_init();
    gmi_sim_start();
    SimUtil_start();
  Sim_readLicenseFile(0);

    Sim_setMessageHandler(messageHandler);
    progress = Progress_new();
    Progress_setCallback(progress, progressHandler);
  }
  catch(pSimError err) {
    printf("Library Initialization Failed\n");
    printf("Did you forget to set SIM_LICENSE_FILE?\n");
    printSimError(err);
    abort();
  }

  //Translate Model
  try {
    // Load the model as an assembly model
    std::string mdlName = argv[1];
    translateModel(mdlName, &simmodel, progress);
  } catch (pSimError err) {
    printSimError(err);
    return(1);
  }

  //Transfer Attributes
  if (argc == 4) {
    try {
      pGModel oldmodel = GM_load(argv[2], NULL, NULL);
      attmngr = SModel_attManager(oldmodel);
      //attmngr = AMAN_load(argv[2]);
      if(!attmngr) {
        fprintf(stderr, "file \"%s\" contains no attributes", argv[2]);
        abort();
      }
      SModel_setAttManager(simmodel, attmngr);
      pACase currentCase;
      pPList allCases = AMAN_cases(attmngr);
      void* iter = 0;
      //Assuming we already have a fairly modern, valid set of 
      //attributes. If converting from old spj, more is needed
      while((currentCase = (pACase) PList_next(allCases, &iter)))
      {
        //update model reference
        printf("Attribute Case Type: %s\n", AttNode_infoType(currentCase));
        printf("Attribute Case Name: %s\n", AttNode_name(currentCase));
        AttCase_setModel(currentCase, simmodel);
      }
      PList_delete(allCases);
      GM_release(oldmodel);
    }
    catch (pSimError err)
    {
      printf("Error Trasferring Attributes\n");
      printSimError(err);
    }
  }

  //Write Output
  try {
    GM_write(simmodel, argv[argc - 1], 0, progress);
  } catch(pSimError err) {
    printf("Problem Outputting\n");
    printSimError(err);
  }

  //Cleanup
  try {
    GM_release(simmodel);
    Progress_delete(progress);

    gmi_sim_stop();
    Sim_unregisterAllKeys();
  SimUtil_stop();
    MS_exit();
  } catch(pSimError err) {
    printf("Cleanup Failed\n");
    printSimError(err);
  }

  return(0);
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
