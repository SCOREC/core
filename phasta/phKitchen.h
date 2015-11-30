#ifndef KITCHEN_H
#define KITCHEN_H

#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <phInput.h>

struct GRStream;
namespace kitchen {
  /** @brief read and attach fields from files */
  void readAndAttachFields(ph::Input& ctrl, apf::Mesh2*& m);
  /** @brief adapt the mesh using the given szFld */
  void adapt(apf::Mesh2* m, apf::Field* szFld);
  /** @brief uniformly refine the mesh */
  void uniformRefinement(ph::Input& ctrl, apf::Mesh2* m);
  /** @brief read fields from the mesh and write to streams */
  void preprocess(gmi_model*& g, apf::Mesh2*& m, ph::Input& ctrl, GRStream* out);
}

#endif
