/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awr.h"
#include "apfMesh.h"
#include "apfField.h"
#include "apfShape.h"

namespace awr {

/* this is pretty simple for now 
   should be robusticized */
void enrichSolution(apf::Field* sol)
{
  apf::Mesh* mesh = sol->getMesh();

  /* solution field information */
  std::string name = apf::getName(sol);
  int type = getValueType(sol);
  apf::FieldShape* shape = sol->getShape();
  int order = shape->getOrder();
  assert (shape == apf::getLagrange(1));
  
  /* enriched solution field information */
  std::string name_e = name + "_e";
  int order_e = order + 1;
  apf::FieldShape* shape_e = apf::getLagrange(order_e);
  apf::Field* sol_e = 
    apf::createField(mesh,name_e.c_str(),type,shape_e);

  /* project old field to new field */
  sol_e->project(sol);
  apf::changeMeshShape(mesh,shape_e,/*project=*/true);

}

}
