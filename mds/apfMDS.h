/******************************************************************************

  Copyright 2014 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef APFMDS_H
#define APFMDS_H

/** \page mds MDS
  The compact Mesh Data Structure is a full-featured parallel unstructured
  mesh representation that is stored mainly in large arrays,
  thus improving locality and reducing memory use.
  It aims to be the default SCOREC mesh data structure.

  A technical description of the storage scheme can be found here:
  http://scorec.rpi.edu/~dibanez/mds.pdf

  Following the design philosophy of APF, we don't even publish the
  internal API for the MDS structure (though it is here and quite clean),
  instead users interact with an MDS mesh through an apf::Mesh2 interface.

  The public API is in apfMDS.h
  */

/** \file apfMDS.h
  \brief Interface to the compact Mesh Data Structure */

#include <map>
//
// AJP: including apf.h for single define: CGNSBCMap
// AJP: alternative is to allow common cgns base header
//      but trying to avoid that since it's not core functionality
#include <apf.h>
//
namespace pcu{
  class PCU;
  PCU* PCU_GetGlobal();
}
struct gmi_model;

namespace apf {

class Mesh;
class Mesh2;
class MeshTag;
class MeshEntity;
class Migration;
class Field;

/** \brief a map from global ids to vertex objects */
typedef std::map<long, MeshEntity*> GlobalToVert;
typedef struct PCU_t PCU_t;


/** \brief create an empty MDS part
  \param model the geometric model interface
  \param dim the eventual mesh dimension. MDS needs to allocate
             arrays based on this before users add entities.
  \param isMatched whether or not there will be matched entities */
Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched, pcu::PCU *PCUObj);

/** \brief load an MDS mesh and model from file
  \param modelfile will be passed to gmi_load to get the model
  \note gmi_register_mesh and gmi_register_null need to be
        called before this function. Also, gmi_register_sim
        may also be called to enable loading of GeomSim, 
        Parasolid, and ACIS models
  */
Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile, pcu::PCU *PCUObj);

/** \brief load an MDS mesh from files
  \param model the geometric model interface
  \param meshfile The path to an SMB format mesh.
                  If the path is "something"".smb", then
                  the file "somethingN.smb" will be loaded
                  where N is the part number.
                  If the path is "something/", then the
                  file "something/N.smb" will be loaded.
                  For both of these cases, if the path is
                  prepended with "bz2:", then it will be uncompressed
                  using PCU file IO functions.
                  Calling apf::Mesh::writeNative on the
                  resulting object will do the same in reverse. */
Mesh2* loadMdsMesh(gmi_model* model, const char* meshfile, pcu::PCU *PCUObj);

// make a serial mesh on all processes - no pmodel & remote link setup
Mesh2* loadSerialMdsMesh(gmi_model* model, const char* meshfile, pcu::PCU *PCUObj);

/** \brief create an MDS mesh from an existing mesh
  \param from the mesh to copy
  \param reorder if true reorder mesh vertices and elements 
         (start from a vertex with minimum Y)
  \param copy_data if true (default), copy Fields/Numberings/Tags
  \details this function uses apf::convert to copy any apf::Mesh */
Mesh2* createMdsMesh(gmi_model* model, Mesh* from, bool reorder=false, bool copy_data=true);

/** \brief apply adjacency-based reordering
  \param t Optional user-defined ordering of the vertices.
           Set this to NULL to use the internal ordering system.
           Otherwise, attach a unique integer to each vertex
           in the range [0, #vertices).
           this will indicate the order in which they appear
           after reordering.
  \details similar to the algorithm for apf::reorder,
           this function will traverse adjacencies to reorder
           each topological type.
           Then all MDS arrays are re-formed in this new order.
           An important side effect of this function is that
           there are no gaps in the MDS arrays after this */
void reorderMdsMesh(Mesh2* mesh, MeshTag* t = 0);

Mesh2* repeatMdsMesh(Mesh2* m, gmi_model* g, Migration* plan, int factor, pcu::PCU *PCUObj);

Mesh2* expandMdsMesh(Mesh2* m, gmi_model* g, int inputPartCount, pcu::PCU *expandedPCU);
/** \brief align the downward adjacencies of matched entities */
bool alignMdsMatches(Mesh2* in);
/** \brief align the downward adjacencies of remote copies */
bool alignMdsRemotes(Mesh2* in);

/** \brief build a null model such that apf::verify accepts the mesh.
  \details given an MDS mesh that is (wrongly) classified on a null model,
  this algorithm will classify all interior entities onto a model region
  and all boundary entities onto a boundary model entity, as defined
  by mesh upward adjacencies. */
void deriveMdsModel(Mesh2* in);

/** \brief Given the mesh vertices that are also model vertices, and the
 *  classification on boundary mesh faces, constructs the classification
 *  on the rest of the boundary entities.
 *
 *  \details Only for tetrahedral mesh with single model region.
 *  The tags provided for face classification are treated as reserved,
 *  and all newly generated tags are distinct regardless of dimension.
 *  It is assumed that the mesh was created using apf::construct, which
 *  by default assigns a tag 0 to the model region. Due to this, it is
 *  advised to provide face tags starting from 1 if uniqueness is desired.
 *  It is assumed that both mesh vertices are indexed from 0 to
 *  (n_verts - 1) and mesh regions from 0 to (n_regions -1).
 *
 *  \param mesh: The mesh in consideration
 *  \param isModelVert Array of bools, one per mesh vertex, telling if that
 *         vertex is also a model vertex
 *  \param nBFaces number of boundary faces
 *  \param bFaces 2D Array of size (n_bfaces x 5). For each face, the row is
 *         [model_face_tag, adj_region_tag, global_vtx_id_1,
 *         global_vtx_id_2, global_vtx_id_3]
 *  \param globalToVert Maps mesh vertex ID to the mesh vertex. Typically
 *         output from apf::construct
 *  \param globalToRegion Maps mesh region ID to the mesh region
 */
void deriveMdlFromManifold(Mesh2* mesh, bool* isModelVert,
			   int nBFaces, int (*bFaces)[5],
			   GlobalToVert &globalToVert,
			   std::map<int, apf::MeshEntity*> &globalToRegion);

/** \brief Given the mesh vertices that are also model vertices, and the
 *  classification on boundary mesh edges, constructs the classification
 *  on the rest of the boundary entities for a 2-D mesh.
 *
 *  \details Only for triangular mesh with single model region.
 *  The tags provided for edge classification are treated as reserved,
 *  and all newly generated tags are distinct regardless of dimension.
 *  It is assumed that both mesh vertices are indexed from 0 to
 *  (n_verts - 1) and mesh faces from 0 to (n_faces -1).
 *
 *  \param mesh: The mesh in consideration
 *  \param isModelVert Array of bools, one per mesh vertex, telling if that
 *         vertex is also a model vertex
 *  \param nBEdges number of boundary faces
 *  \param bEdges 2D Array of size (nBEdges x 4). For each face, the row is
 *         [model_edge_tag, adj_face_tag, global_vtx_id_1, global_vtx_id_2]
 *  \param globalToVert Maps mesh vertex ID to the mesh vertex. Typically
 *         output from apf::construct
 *  \param globalToFace Maps mesh face ID to the mesh face
 */
void derive2DMdlFromManifold(Mesh2* mesh, bool* isModelVert,
			     int nBEdges, int (*bEdges)[4],
			     GlobalToVert &globalToVert,
			     std::map<int, apf::MeshEntity*> &globalToFace);

/** \brief change the dimension of an MDS mesh
  \details this should be called before adding entities of
  dimension higher than the previous mesh dimension
  (when building a higher dimensional mesh from a lower one),
  or after removing all entities of higher dimension
  (when reducing a high dimensional mesh to a lower one) */
void changeMdsDimension(Mesh2* in, int d);

/** \brief returns the dimension-unique index for this entity
 \details this function only works when the arrays have no gaps,
 so call apf::reorderMdsMesh after any mesh modification. */
int getMdsIndex(Mesh2* in, MeshEntity* e);

/** \brief retrieve an entity by dimension and index
  \details indices follow iteration order, so this
  function is equivalent to iterating (index) times,
  but is actually much faster than that.
  this function only works when the arrays have no gaps,
  so call apf::reorderMdsMesh after any mesh modification. */
MeshEntity* getMdsEntity(Mesh2* in, int dimension, int index);

Mesh2* loadMdsFromCGNS(gmi_model* g, const char* filename, CGNSBCMap& cgnsBCMap);
Mesh2* loadMdsFromCGNS2(PCU_t h, gmi_model* g, const char* filename, CGNSBCMap& cgnsBCMap);

// names of mesh data to read from file: (VERTEX, VelocityX; CellCentre, Pressure)
Mesh2* loadMdsFromCGNS(gmi_model* g, const char* filename, CGNSBCMap& cgnsBCMap, const std::vector<std::pair<std::string, std::string>>& meshData);
Mesh2* loadMdsFromCGNS2(PCU_t h, gmi_model* g, const char* filename, CGNSBCMap& cgnsBCMap, const std::vector<std::pair<std::string, std::string>>& meshData);


int gmshMajorVersion(const char* filename);

Mesh2* loadMdsFromGmsh(gmi_model* g, const char* filename, pcu::PCU *PCUObj);

Mesh2* loadMdsDmgFromGmsh(const char* fnameDmg, const char* filename, pcu::PCU *PCUObj);

Mesh2* loadMdsFromUgrid(gmi_model* g, const char* filename, pcu::PCU *PCUObj);

void printUgridPtnStats(gmi_model* g, const char* ugridfile, const char* ptnfile,
    const double elmWeights[], pcu::PCU *PCUObj);

/** \brief load an MDS mesh from ANSYS .node and .elem files
  \details this call takes two filenames, one
  for a .node and another for a .elem file.

  the resulting MDS mesh will be constructed with a null
  geometric model via gmi_load(".null"), so be sure to
  call gmi_register_null before this function.

  currently, ANSYS element types SOLID72 and SOLID92 are
  supported, which become linear and quadratic tetrahedra,
  respectively. */
Mesh2* loadMdsFromANSYS(const char* nodefile, const char* elemfile, pcu::PCU *PCUObj);

void disownMdsModel(Mesh2* in);

void setMdsMatching(Mesh2* in, bool has);

Mesh2* loadMdsPart(gmi_model* model, const char* meshfile, pcu::PCU *PCUObj);
void writeMdsPart(Mesh2* m, const char* meshfile);

}

#endif
