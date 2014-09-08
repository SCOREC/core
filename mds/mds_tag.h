/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MDS_TAG_H
#define MDS_TAG_H

#include "mds.h"

struct mds_tag {
  struct mds_tag* next;
  int bytes;
  int user_type;
  char* data[MDS_TYPES];
  unsigned char* has[MDS_TYPES];
  char* name;
};

struct mds_tags {
  struct mds_tag* first;
};

void mds_create_tags(struct mds_tags* ts);
void mds_destroy_tags(struct mds_tags* ts);
void mds_grow_tags(
    struct mds_tags* ts,
    struct mds* m,
    mds_id old_cap[MDS_TYPES]);
struct mds_tag* mds_create_tag(
    struct mds_tags* ts,
    const char* name,
    int bytes,
    int user_type);
void mds_destroy_tag(struct mds_tags* ts, struct mds_tag* t);
void* mds_get_tag(struct mds_tag* tag, mds_id e);
struct mds_tag* mds_find_tag(struct mds_tags* ts, const char* name);
int mds_has_tag(struct mds_tag* tag, mds_id e);
void mds_give_tag(struct mds_tag* tag, struct mds* m, mds_id e);
void mds_take_tag(struct mds_tag* tag, mds_id e);

void mds_swap_tag_structs(struct mds_tags* as, struct mds_tag** a,
    struct mds_tags* bs, struct mds_tag** b);

#endif
