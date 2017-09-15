/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef _PUMI_ERROR_CODE_H_
#define _PUMI_ERROR_CODE_H_

#ifdef __cplusplus
extern "C" {

#else   /* not __cplusplus */

#endif /* not __cplusplus */

 enum PUMI_ErrCode { /* general or common errors */
  PUMI_SUCCESS = 0, // no error
  PUMI_MESH_ALREADY_LOADED,
  PUMI_FILE_NOT_FOUND,
  PUMI_FILE_WRITE_ERROR,
  PUMI_NULL_ARRAY,
  PUMI_BAD_ARRAY_SIZE,
  PUMI_BAD_ARRAY_DIMENSION,
  PUMI_INVALID_ENTITY_HANDLE,
  PUMI_INVALID_ENTITY_COUNT,
  PUMI_INVALID_ENTITY_TYPE,
  PUMI_INVALID_ENTITY_TOPO,
  PUMI_BAD_TYPE_AND_TOPO,
  PUMI_ENTITY_CREATION_ERROR,
  PUMI_INVALID_TAG_HANDLE,
  PUMI_TAG_NOT_FOUND,
  PUMI_TAG_ALREADY_EXISTS,
  PUMI_TAG_IN_USE, // try to delete a tag that is in use
  PUMI_INVALID_SET_HANDLE,
  PUMI_INVALID_ITERATOR,
  PUMI_INVALID_ARGUMENT,
  PUMI_MEMORY_ALLOCATION_FAILED,
  PUMI_NOT_SUPPORTED,
  PUMI_FAILURE, 
  // scorec-defined error codes from this line, shouldn't be used in itaps impl
  PUMI_INVALID_MESH_INSTANCE,
  PUMI_INVALID_GEOM_MODEL,
  PUMI_INVALID_GEOM_CLAS,
  PUMI_INVALID_PTN_CLAS,
  PUMI_INVALID_REMOTE,
  PUMI_INVALID_MATCH,
  PUMI_INVALID_PART_HANDLE, 
  PUMI_INVALID_PART_ID,
  PUMI_INVALID_SET_TYPE,
  PUMI_INVALID_TAG_TYPE,
  PUMI_ENTITY_NOT_FOUND,
  PUMI_ENTITY_ALREADY_EXISTS,
  PUMI_REMOTE_NOT_FOUND,
  PUMI_GHOST_NOT_FOUND,
  PUMI_CB_ERROR
  };

#ifdef __cplusplus
}
#endif

#endif  /* ifndef _PUMI_H_ */


