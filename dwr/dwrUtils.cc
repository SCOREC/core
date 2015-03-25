/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <cstdarg>
#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include "dwrUtils.h"

namespace dwr {

static double getEdgeLength(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::MeshElement* element = apf::createMeshElement(m,e);
  double h = apf::measure(element);
  apf::destroyMeshElement(element);
  return h;
}

double getMeshSize(apf::Mesh* m, apf::MeshEntity* e)
{
  /* right now the maximum edge length */
  double h = 0.0;
  apf::Downward edges;
  int ne = m->getDownward(e,1,edges);
  for (int i=0; i < ne; ++i)
    h = std::max(h, getEdgeLength(m,edges[i]));
  return h;
}

void print(const char* format, ...)
{
  if (PCU_Comm_Self())
    return;
  printf("DWR: ");
  va_list ap;
  va_start(ap,format);
  vfprintf(stdout,format,ap);
  va_end(ap);
  printf("\n");
}

}
