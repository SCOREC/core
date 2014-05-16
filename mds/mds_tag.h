/*
   Copyright 2014 Dan Ibanez

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

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
    struct mds* m,
    const char* name,
    int bytes,
    int user_type);
void mds_destroy_tag(struct mds_tags* ts, struct mds_tag* t);
void* mds_get_tag(struct mds_tag* tag, mds_id e);
struct mds_tag* mds_find_tag(struct mds_tags* ts, const char* name);
int mds_has_tag(struct mds_tag* tag, mds_id e);
void mds_give_tag(struct mds_tag* tag, struct mds* m, mds_id e);
void mds_take_tag(struct mds_tag* tag, mds_id e);

#endif
