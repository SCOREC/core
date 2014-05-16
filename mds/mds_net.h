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

#ifndef MDS_NET_H
#define MDS_NET_H

#include "mds.h"

struct mds_copy {
  mds_id e;
  int p;
};

struct mds_copies {
  int n;
  struct mds_copy c[];
};

struct mds_net {
  mds_id n[MDS_TYPES];
  struct mds_copies** data[MDS_TYPES];
};

struct mds_links {
  unsigned np;
  unsigned* n;
  unsigned* p;
  unsigned** l;
};

void mds_create_net(struct mds_net* net);
void mds_destroy_net(struct mds_net* net, struct mds* m);
struct mds_copies* mds_make_copies(int n);
void mds_set_copies(struct mds_net* net, struct mds* m, mds_id e,
    struct mds_copies* c);
struct mds_copies* mds_get_copies(struct mds_net* net, mds_id e);
void mds_grow_net(
    struct mds_net* net,
    struct mds* m,
    mds_id old_cap[MDS_TYPES]);

void mds_add_copy(struct mds_net* net, struct mds* m, mds_id e,
    struct mds_copy c);

void mds_get_type_links(struct mds_net* net, struct mds* m,
    int t, struct mds_links* ln);
void mds_set_type_links(struct mds_net* net, struct mds* m,
    int t, struct mds_links* ln);
void mds_free_links(struct mds_links* ln);

int mds_net_empty(struct mds_net* net);

void mds_get_local_matches(struct mds_net* net, struct mds* m,
                         int t, struct mds_links* ln);
void mds_set_local_matches(struct mds_net* net, struct mds* m,
                         int t, struct mds_links* ln);
void mds_free_local_links(struct mds_links* ln);

#endif
