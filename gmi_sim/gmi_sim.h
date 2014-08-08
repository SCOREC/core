/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_SIM_H
#define GMI_SIM_H

class SGModel;

void gmi_register_sim(void);
struct gmi_model* gmi_import_sim(SGModel* sim_model);
SGModel* gmi_export_sim(struct gmi_model* m);

#endif


