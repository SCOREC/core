/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_CAP_H
#define GMI_CAP_H


#include "gmi.h"
#ifdef __cplusplus
#include <string>
#include <vector>
#include "CreateMG_Framework_Geometry.h"
#else
extern "C" {
#endif

void gmi_cap_start(void);
void gmi_cap_stop(void);
void gmi_register_cap(void);
struct gmi_model* gmi_cap_load(const char* creFileName);

#ifdef __cplusplus

void gmi_cap_probe(
  const char* creFileName, std::vector<std::string>& mesh_names
);
void gmi_cap_probe(
  const char* creFileName, std::string& model_content,
  std::vector<std::string>& mesh_names, std::vector<std::string>& mesh_contents
);
struct gmi_model* gmi_cap_load_some(
  const char* creFileName, const std::vector<std::string>& mesh_names
);

gmi_ent* toGmiEntity(CreateMG::M_GTopo topo);
CreateMG::M_GTopo fromGmiEntity(gmi_ent* g);

struct gmi_model* gmi_import_cap(CreateMG::GDBI* gi);
struct gmi_model* gmi_import_cap(CreateMG::GDBI* gi, CreateMG::M_GModel);
CreateMG::GDBI* gmi_export_cap(struct gmi_model* m);

#endif // __cplusplus


#ifndef __cplusplus
} // extern "C"
#endif

#endif // GMI_CAP_H

