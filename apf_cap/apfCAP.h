/****************************************************************************** 

  Copyright 2025 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

******************************************************************************/
#ifndef APF_CAP_H
#define APF_CAP_H

#include <string>

// Forward declarations
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

namespace apf {

class Mesh2;
class Field;

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
 * gmi_cap_load, or gmi_cap_load_some. Try to load the first mesh (by index).
 *
 * \param model A gmi_model associated with a Capstone geometry.
 * \return an apf::Mesh2 interface to the Capstone mesh.
 */
Mesh2* createCapMesh(gmi_model* model, pcu::PCU* PCUObj);

/**
 * \brief Create an apf::Mesh2 object from a mesh linked to a loaded
 * geometry model.
 *
 * The gmi_model should be loaded previously by gmi_load (on a .cre file),
 * gmi_cap_load, or gmi_cap_load_some. The list of acceptable values for
 * meshname can be found by using gmi_cap_probe or direct Capstone interfaces.
 *
 * \param model A gmi_model associated with a Capstone geometry.
 * \param meshname The name of a mesh associated with
 * \return an apf::Mesh2 interface to the Capstone mesh.
 */
Mesh2* createCapMesh(gmi_model* model, const char* meshname, pcu::PCU* PCUObj);

/**
 * \brief Generate Capstone mesh object on a model linked to a Capstone model.
 *
 * \param model A Capstone GMI model to generate a mesh on.
 * \param meshname Name for the new generated mesh.
 * \param PCUObj The PCU object to link to the mesh.
 * \return An apf::Mesh2 object with the mesh.
 * \bug dimension != 3 is not supported.
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

}//namespace apf

#endif // APF_CAP_H
