#ifndef CHEF_H
#define CHEF_H

/** \file chef.h
    \brief The Chef interface */

#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <phInput.h>

struct RStream;
struct GRStream;
namespace chef {
  /** @brief read and write to and from files */
  void cook(gmi_model*& g, apf::Mesh2*& m, pcu::PCU *PCUObj);
  /** @brief read and write to and from files and load control info*/
  void cook(gmi_model*& g, apf::Mesh2*& m, ph::Input& ctrl, pcu::PCU *PCUObj);
  /** @brief read from stream and write to files */
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* in, pcu::PCU *PCUObj);
  /** @brief read from files and write to stream */
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, GRStream* out, pcu::PCU *PCUObj);
  /** @brief read and write to and from streams */
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* in, GRStream* out, pcu::PCU *PCUObj);
  /** @brief extract a field from a packed field */
  apf::Field* extractField(apf::Mesh* m,
    const char* packedFieldname,
    const char* requestFieldname,
    int firstComp,
    int numOfComp,
    bool simField = false);
  /** @brief combine 3 fields into 1 packed field */
  apf::Field* combineField(apf::Mesh* m,
    const char* packedFieldname,
    const char* inFieldname1,
    const char* inFieldname2,
    const char* inFieldname3);
  /** @brief read and attach fields from files */
  void readAndAttachFields(ph::Input& ctrl, apf::Mesh2*& m);
  /** @brief load balance the partition then reorder the vertices */
  void balanceAndReorder(ph::Input& ctrl, apf::Mesh2* m);
  /** @brief load balance the partition */
  void balance(ph::Input& ctrl, apf::Mesh2* m);
  /** @brief adapt the mesh using the given szFld */
  void adapt(apf::Mesh2* m, apf::Field* szFld);
  /** @brief adapt the mesh based on input control */
  void adapt(apf::Mesh2* m, apf::Field* szFld, ph::Input& ctrl);
  /** @brief uniformly refine the mesh */
  void uniformRefinement(ph::Input& ctrl, apf::Mesh2* m);
  /** @brief adapt around the zero level set */
  void adaptLevelSet(ph::Input& ctrl, apf::Mesh2* m);
  /** @brief read fields from the mesh and write to files */
  void preprocess(apf::Mesh2*& m, ph::Input& ctrl);
  /** @brief read fields from the mesh and write to streams */
  void preprocess(apf::Mesh2*& m, ph::Input& ctrl, GRStream* out);
}

#endif
