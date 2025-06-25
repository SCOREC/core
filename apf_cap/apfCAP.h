/****************************************************************************** 

  Copyright 2025 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

******************************************************************************/
#ifndef APF_CAP_H
#define APF_CAP_H
/**
 * \file apfCAP.h
 * \brief Capstone apf::Mesh2 implementation and interface.
 *
 * Like gmi_cap, the interface is used in two ways:
 *
 * 1. to import an existing Capstone mesh; and
 * 2. to load or generate a mesh associated with a model which was previously
 *    loaded by gmi_cap.
 *
 * \note Files which `#include` apfCAP.h should also `#include`
 * CreateMG_Framework_Mesh.h (or another Capstone header with the full
 * MeshDatabaseInterface definition) to import meshes or use getCapNative.
 */

#include <string>

/**
 * \cond
 * Forward declarations
 */
struct gmi_model;
namespace pcu {
  class PCU;
}
namespace CreateMG {
  namespace Geometry { class GeometryDatabaseInterface; }
  namespace Mesh { class MeshDatabaseInterface; }
  typedef Geometry::GeometryDatabaseInterface GDBI;
  typedef Mesh::MeshDatabaseInterface MDBI;
}
/** \endcond */

namespace apf {

class Mesh2;
class Field;

/**
 * \defgroup apf_cap Capstone APF mesh interface
 *
 * apf_cap provides access to an implementation of apf::Mesh2. The model must
 * be loaded by gmi_cap. The interface in apfCAP.h provides additional
 * functions to simplify loading and interacting with the underlying mesh.
 *
 * \{
 */

/**
 * \brief Test for compiled Capstone library support.
 * \return true if apf_cap was compiled with Capstone. otherwise all routines
 *         will call apf::fail.
 */
bool hasCAP() noexcept;

/**
 * \brief Create an apf::Mesh2 object from a Capstone mesh database.
 *
 * This object should be destroyed by apf::destroyMesh. Since Capstone meshes
 * are serial only right now, PCUObj will not be directly used, but may have an
 * impact on collective calls involving the mesh (e.g. a PCU::Add which uses
 * Mesh::getPCU()).
 *
 * \param mdb A Capstone MeshDatabaseInterface with a current M_MModel.
 * \param gdb A Capstone GeometryDatabaseInterface with a current M_GModel.
 * \param PCUObj The PCU communicator to define the mesh over.
 */
Mesh2* createCapMesh(
  CreateMG::MDBI* mdb, CreateMG::GDBI* gdb, pcu::PCU *PCUObj
);

/**
 * \brief Create an apf::Mesh2 object from the first mesh linked to a loaded
 * geometry model.
 *
 * The gmi_model should be loaded previously by gmi_load (on a .cre file),
 * gmi_cap_load, or gmi_cap_load_selective. Try to load the first mesh (by
 * index).
 *
 * \param model A gmi_model associated with a Capstone geometry.
 * \param PCUObj The PCU communicator to define the mesh over
 * \return an apf::Mesh2 interface to the Capstone mesh.
 */
Mesh2* createCapMesh(gmi_model* model, pcu::PCU* PCUObj);

/**
 * \brief Create an apf::Mesh2 object from a mesh linked to a loaded
 * geometry model.
 *
 * The gmi_model should be loaded previously by gmi_load (on a .cre file),
 * gmi_cap_load, or gmi_cap_load_selective. The list of acceptable values for
 * meshname can be found by using gmi_cap_probe or direct Capstone interfaces.
 *
 * \param model A gmi_model associated with a Capstone geometry.
 * \param meshname The name of a mesh associated with
 * \param PCUObj The PCU communicator to define the mesh over
 * \return an apf::Mesh2 interface to the Capstone mesh.
 */
Mesh2* createCapMesh(gmi_model* model, const char* meshname, pcu::PCU* PCUObj);

/**
 * \brief Generate Capstone mesh object on a model linked to a Capstone model.
 *
 * \param model A Capstone GMI model to generate a mesh on.
 * \param dimension The dimension of mesh to generate
 * \param PCUObj The PCU object to link to the mesh.
 * \return An apf::Mesh2 object with the mesh.
 */
Mesh2* generateCapMesh(
  gmi_model* model, int dimension, pcu::PCU* PCUObj
);

/**
 * \brief Make an empty Capstone M_MModel associated with the model.
 *
 * \param model Previously loaded Capstone gmi_model.
 * \param meshname Name for new mesh in the Capstone database.
 * \param PCUObj The PCU object to associate the new apf::Mesh2 with.
 * \return a new apf::Mesh2 object.
 */
Mesh2* makeEmptyCapMesh(
  gmi_model* model, const char* meshname, pcu::PCU* PCUObj
);

/**
 * Disown capMesh's gmi_model.
 *
 * Mark the gmi_model as non-owned so that the destructor does not call
 * gmi_destroy.
 *
 * \param capMesh A Capstone mesh wrapper.
 */
void disownCapModel(Mesh2* capMesh);

/**
 * \brief Get native Capstone mesh database interface.
 * \param capMesh Previously loaded capstone mesh.
 * \return Underlying Capstone mesh database interface from capMesh.
 */
CreateMG::MDBI* getCapNative(Mesh2* capMesh);

/**
 * \brief Test for smoothCAPAnisoSizes support.
 *
 * \return A boolean indicating whether support was compiled. False indicates
 *         the call would fail.
 *
 * \details smoothCAPAnisoSizes is only compiled if support for the underlying
 *          call is detected in the version of Capstone apf_cap was compiled
 *          against. Otherwise the call will always apf::fail. Use this
 *          function to programmatically test for the capability.
 */
bool has_smoothCAPAnisoSizes(void) noexcept;

/**
 * \brief Use the SizingMetricTool to smooth a size field on a Capstone mesh.
 *
 * \param m A Capstone mesh.
 * \param analysis The Capstone analysis to use.
 * \param frames An apf::Field of apf::Matrix3x3 with orthogonal basis frames.
 * \param scales An apf::Field of apf::Vector3 with frame scales (eigenvalues).
 * \return A boolean indicating success.
 * \pre m must be a Capstone mesh.
 */
bool smoothCAPAnisoSizes(
  apf::Mesh2* m, std::string analysis, apf::Field* scales, apf::Field* frames
);

/** \} */

} // namespace apf

#endif // APF_CAP_H
