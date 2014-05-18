/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_SIM_H
#define GMI_SIM_H

class SGModel;

extern "C" void gmi_register_sim(void);
extern "C" struct gmi_model* gmi_import_sim(SGModel* sim_model);

#endif


