/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MDS_NET_H
#define MDS_NET_H

#include "mds.h"

struct mds_copy {
  mds_id e;
  int p;
};

struct mds_copies {
  int n;
  struct mds_copy c[1];
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
#define MDS_LINKS_INIT {0,0,0,0}

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
