#include <PCU.h>
#include <lionPrint.h>
#include "phGrowthCurves.h"
#include "phOutput.h"
namespace ph {
  void getGrowthCurves(Output& o) {
    o.nGrowthCurves = 0;
    o.nLayeredMeshVertices = 0;
    if(PCU_Comm_Self() == 0)
      lion_oprint(1,"warning! \'%s\' requires the Simmetrix SimAdvMeshing library\n",__func__);
    return;
  }
}
