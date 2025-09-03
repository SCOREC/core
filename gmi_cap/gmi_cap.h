/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_CAP_H
#define GMI_CAP_H
/**
 * \file gmi_cap.h
 * \brief Concrete Capstone GMI implementation.
 *
 * This interface can be used in two ways:
 *
 * 1. to import an existing Capstone model; and
 * 2. to load the model and (optionally) meshes from a CRE file by wrapping
 *    Capstone calls.
 *
 * The first option results in a non-owned model (gmi_cap_start does not need
 * to be called). The second option results in an owned model (gmi_cap_start
 * MUST be called) for which destruction is handled by gmi_destroy. Both
 * options require gmi_register_cap to be called. It is recommended to call
 * gmi_cap_start (if desired) before gmi_register_cap.
 */

#ifdef __cplusplus
#include <string>
#include <vector>
#include <CreateMG_CoreDefs.h>
// Forward declaration
namespace CreateMG {
  namespace Geometry { class GeometryDatabaseInterface; }
  typedef Geometry::GeometryDatabaseInterface GDBI;
}
#else
extern "C" {
#endif
// Forward declarations
struct gmi_model;
struct gmi_ent;

/**
 * \defgroup gmi_cap Capstone geometric model interface
 *
 * gmi_cap implements gmi. The interface in gmi_cap.h includes additional
 * functions to initialize the library and simplify loading of Capstone models.
 *
 * \{
 */

/**
 * \brief Initialize gmi_cap library.
 *
 * \note This call is required before calling any other gmi_cap functions.
 */
void gmi_cap_start(void);

/**
 * \brief Returns true if gmi_cap library has been initialized with gmi_cap_start
 */
bool is_gmi_cap_started();

/**
 * \brief Finalize gmi_cap library.
 *
 * \note No gmi_cap functions except gmi_cap_import and gmi_cap_export should
 * bec called after this function.
 */
void gmi_cap_stop(void);
/**
 * \brief Register gmi_cap in GMI's global list.
 *
 * \note This call is required before every gmi_cap function except
 * gmi_cap_start.
 *
 * This call registers gmi_cap with GMI. It also registers the ".cre" extension
 * with the Create loader so that gmi_load can be used with ".cre" files.
 */
void gmi_register_cap(void);
/**
 * \brief Load a gmi_model from a CRE file.
 *
 * This method is called by gmi_load when the input file has a ".cre"
 * extension. It loads the file with all associated mesh models. Use
 * gmi_cap_load_selective to only load certain meshes.
 *
 * \param creFileName the name of the CRE file
 * \pre gmi_cap_stop must be called before this method.
 * \pre creFileName ends in ".cre"
 */
struct gmi_model* gmi_cap_load(const char* creFileName);
/**
 * \brief Write a Capstone gmi_model and associated mesh to a CRE file.
 *
 * \param model A previously loaded gmi_model from gmi_cap.
 * \param creFileName filename for the new CRE file.
 */
void gmi_cap_write(struct gmi_model* model, const char* creFileName);

#ifdef __cplusplus

/**
 * \brief Probe a CRE file for associated mesh names.
 *
 * This function calls Capstone APIs to probe the CRE file header for meshes.
 * It is intended to be used in conjunction with gmi_cap_load_selective.
 *
 * \param[in]  creFileName CRE file to probe
 * \param[out] mesh_names A vector to be filled with mesh names
 */
void gmi_cap_probe(
  const char* creFileName, std::vector<std::string>& mesh_names
);
/**
 * \brief Probe a CRE file for associated mesh names and contents.
 *
 * This function calls Capstone APIs to probe the CRE file header for model
 * contents, mesh names, and mesh contents. It is intended to be used in
 * conjunction with gmi_cap_load_selective, and/or to assist when there are
 * multiple loadable meshes (i.e. by examining contents or displaying an
 * interactive menu).
 *
 * \warning The format of model_content and mesh_contents is governed by the
 * Capstone library and we cannot make any claims about its stability.
 *
 * \param[in]  creFileName CRE file to probe
 * \param[out] model_content A string indicating mesh contents
 * \param[out] mesh_names A vector to be filled with mesh names
 * \param[out] mesh_contents A vector to be filled with strings indicating mesh
 *                           contents
 */
void gmi_cap_probe(
  const char* creFileName, std::string& model_content,
  std::vector<std::string>& mesh_names, std::vector<std::string>& mesh_contents
);
/**
 * \brief Load a gmi_model from a CRE file and select the associated loaded
 * meshes.
 *
 * Unlike other data formats, CRE files contain both geometries and meshes. As
 * such, the meshes must be selected at load time. The normal gmi_cap_load
 * (which is registered to gmi_load with the ".cre" extension) loads all
 * meshes by default. In certain scenarios this is undesirable (e.g. when you
 * only want the first mesh in a CRE file multiple meshes, or when you only
 * want to load the mesh on one rank but need the geometry on all ranks).
 *
 * \param creFileName The CRE file to load
 * \param mesh_names A vector of mesh names to load, which may be empty.
 * \return A new gmi_model linked to a Capstone geometry model
 */
struct gmi_model* gmi_cap_load_selective(
  const char* creFileName, const std::vector<std::string>& mesh_names
);

gmi_ent* toGmiEntity(CreateMG::M_GTopo topo);
CreateMG::M_GTopo fromGmiEntity(gmi_ent* g);

/**
 * \brief Import a non-owned Capstone model.
 *
 * Loads the current geometry model on the interface.
 *
 * \param gi Geometry interface to load.
 * \pre gi has a current model.
 */
struct gmi_model* gmi_import_cap(CreateMG::GDBI* gi);
/**
 * \brief Import a non-owned, specific Capstone model.
 *
 * \param gi Geometry interface to load from.
 * \param gmodel The existing geometry model to load.
 */
struct gmi_model* gmi_import_cap(
  CreateMG::GDBI* gi, CreateMG::M_GModel gmodel
);
/**
 * \brief Expose the underlying Capstone geometry interface from a gmi_model.
 *
 * \param m Loaded or imported gmi_model.
 */
CreateMG::GDBI* gmi_export_cap(struct gmi_model* m);

#endif // __cplusplus


#ifndef __cplusplus
} // extern "C"
#endif

/** \} */

#endif // GMI_CAP_H

