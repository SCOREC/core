#include "MeshSim.h"
#include "SimAttribute.h"
#include "AttributeTypes.h"
#include "SimParasolidKrnl.h"
#include "string.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

/* Secret Undocumented Simmetrix Functions! */
void SModel_setAttManager(pModel model, pAManager attMan);

int main(int argc, char* argv[])
{
  if (argc != 4) {
    printf("Usage: %s <Parasolid model> <Simmetrix .spj> <output .smd>\n", argv[0]);
    printf("       to convert the Simmetrix .spj file into a .smd file\n");
    return 0;
  }

  const char* gname = argv[1];
  const char* fname = argv[2];
  const char* oname = argv[3];
  pGModel model;

  MS_init();
  SimModel_start();
  SimParasolid_start(1);
  SimUtil_start();
  Sim_readLicenseFile(0);

  pProgress prog = Progress_new();
  Progress_setDefaultCallback(prog);

  // create Parasolid native model from part file
  pNativeModel pnModel;
  pnModel = ParasolidNM_createFromFile(gname,0); 
  if(NM_isAssemblyModel(pnModel)) {
    pGAModel amodel = GAM_createFromNativeModel(pnModel, prog);
    NM_release(pnModel);
    model = GM_createFromAssemblyModel(amodel, NULL, prog);
    GM_release(amodel);
    pnModel = GM_nativeModel(model);
  }
  else
    model = GM_createFromNativeModel(pnModel, prog);

  // initialize the attribute manager
  pAManager attmngr = AMAN_load(fname);
  if (attmngr == 0){
    std::cout << "Error: could not open attribute file "<<fname << std::endl;
    exit(-1);
  }

  // Associate the attribute manager with the model and write out both the native and the smd
  SModel_setAttManager(model,attmngr);
  // Here we need to fix up the image for the cases to be "phasta"
  // to match the new attDefs, and set the model
  pPList allCases = AMAN_cases(attmngr);


  pACase currentCase;
  void * iter = 0;
  while ((currentCase = (pACase)PList_next(allCases,&iter))) {
    char* caseInfoType = AttNode_infoType(currentCase);
    char* caseName = AttNode_name(currentCase);
    std::cout << "Info Type for this case: "<<caseInfoType << std::endl;
    std::cout << "Name for this case: "<<caseName << std::endl;
    // delete the meshing case
    if(!strcmp(caseInfoType,"meshing")) {
      // The meshing case doesn't seem to get fully deleted, so still want
      // to set the model so the smd file has the correct path
      AttCase_setModel(currentCase,model);
      std::cout << "Deleting case "<<caseInfoType << std::endl;;
      AMAN_removeNode(attmngr,currentCase,1);
    } else {
      std::cout << "Current image for class is "<<AttNode_imageClass(currentCase) << std::endl;
      AttNode_setImageClass(currentCase,"phasta");
      std::cout << "Now set to "<<AttNode_imageClass(currentCase) << std::endl;
      std::cout << "Now set the model\n" << std::endl;
      AttCase_setModel(currentCase,model);

      //see if analysis case exists
      if(!AMAN_findCaseByType(attmngr, "analysis")) {
        printf("no analysis case found, creating one\n");
        pACase analysis = AMAN_newCase(attmngr, "analysis","analysis",(pModel)model);
        AttNode_setImageClass(analysis,"phasta");
        AttNode_add(analysis, currentCase);
        AttCase_setModel(analysis,model);
        AttCase_setVersion(analysis, 1);
        printf("set att. version to 1 for analysis case, added child\n");
      }
      // For analysis, also set the version of the attdefs. 
      // This is the number in the *.adi file in []
      if(!strcmp(caseInfoType,"analysis")) {
        std::cout << "Setting version of attributes for case "<<caseInfoType<<" to 1" << std::endl;
        AttCase_setVersion(currentCase,1);
      }
    }
    Sim_deleteString(caseInfoType);
    Sim_deleteString(caseName);
  }
  PList_delete(allCases);

  // Now write out the files, these can be read into SimModeler
  GM_write(model, oname, 0, prog);

  NM_release(pnModel);
  GM_release(model);

  AMAN_release(attmngr);

  SimModel_stop();
  SimParasolid_stop(1);
  MS_exit();
}   
