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

struct gmi_model;

namespace apf {

class Mesh;
class Mesh2;
class MeshTag;
class MeshEntity;
class Migration;

/** \brief create an empty MDS part
  \param model the geometric model interface
  \param dim the eventual mesh dimension. MDS needs to allocate
             arrays based on this before users add entities.
  \param isMatched whether or not there will be matched entities */
Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched);

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
Mesh2* loadMdsMesh(gmi_model* model, const char* meshfile);

/** \brief load an MDS mesh and model from file
  \param modelfile will be passed to gmi_load to get the model */
Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile);

/** \brief create an MDS mesh from an existing mesh
  \param from the mesh to copy
  \details this function uses apf::convert to copy any apf::Mesh */
Mesh2* createMdsMesh(gmi_model* model, Mesh* from);

/** \brief apply adjacency-based reordering
  \details similar to the algorithm for apf::reorder,
           this function will traverse adjacencies to reorder
           each topological type.
           Then all MDS arrays are re-formed in this new order.
           An important side effect of this function is that
           there are no gaps in the MDS arrays after this */
void reorderMdsMesh(Mesh2* mesh);

/** \brief split an MDS mesh into multiple parts per process using threads
  \param m the MDS mesh
  \param plan the plan for splitting the local part, usually the output
              of some apf::Splitter
  \param n the number of threads per process. this must be consistent
           on all processes.
  \param runAfter the code to run on the split parts. This function is
                  called by n threads, each with a unique MDS split part.
  */
void splitMdsMesh(Mesh2* m, Migration* plan, int n, void (*runAfter)(Mesh2*));
Mesh2* repeatMdsMesh(Mesh2* m, gmi_model* g, Migration* plan, int factor);

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

Mesh2* loadMdsFromGmsh(gmi_model* g, const char* filename);

/** \brief load an MDS mesh from ANSYS .node and .elem files
  \details this call takes two filenames, one
  for a .node and another for a .elem file.

  the resulting MDS mesh will be constructed with a null
  geometric model via gmi_load(".null"), so be sure to
  call gmi_register_null before this function.

  currently, ANSYS element types SOLID72 and SOLID92 are
  supported, which become linear and quadratic tetrahedra,
  respectively. */
Mesh2* loadMdsFromANSYS(const char* nodefile, const char* elemfile);

}

#endif
