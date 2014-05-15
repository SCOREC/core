/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_COLL_H
#define PCU_COLL_H

#include "pcu_mpi.h"

/* The PCU Collectives system (pcu_coll for short) implements
   non-blocking collective operations based loosely on binary-tree
   or binomial communication patterns.

   This system is an abstraction and implementation of reduction, broadcast,
   and scan algorithms.
   Because all communication uses the pcu_mpi primitives, the
   system works in hybrid mode as well. */

/* The pcu_merge is the equivalent of the MPI_Op.
   arguments are usually arrays of some type,
   and the operations is sum, min, max, etc. */

typedef void pcu_merge(void* local, void* incoming, size_t size);
void pcu_merge_assign(void* local, void* incoming, size_t size);
void pcu_add_doubles(void* local, void* incoming, size_t size);
void pcu_max_doubles(void* local, void* incoming, size_t size);
void pcu_min_doubles(void* local, void* incoming, size_t size);
void pcu_add_ints(void* local, void* incoming, size_t size);
void pcu_min_ints(void* local, void* incoming, size_t size);
void pcu_max_ints(void* local, void* incoming, size_t size);
void pcu_add_longs(void* local, void* incoming, size_t size);

/* Enumerated actions that a rank takes during one
   step of the communication pattern */
enum
{
  pcu_coll_send,
  pcu_coll_recv,
  pcu_coll_idle
};

/* The pcu_pattern is an abstraction of a communication
   pattern that takes O(lg(n)) steps for n peers,
   and at each step (rank) communicates with (rank +- 2^k)
   where k is some integer.
   The entire pattern state is encoded into this integer k, which
   for convenience is stored as bit=2^k=1<<k.
 */
typedef struct
{
  int (*begin_bit)(void); //initialize state bit
  bool (*end_bit)(int bit); //return true if bit is one past the last
  int (*action)(int bit); //return action enum for this step
  int (*peer)(int bit); //return the peer to communicate with
  int (*shift)(int bit); //shift the bit up or down
} pcu_pattern;

/* The pcu_coll object stores the state of a non-blocking
   collective operation. */
typedef struct
{
  pcu_pattern* pattern; //communication pattern controller
  pcu_merge* merge; //merge operation
  pcu_message message; //local data being operated on
  int bit; //pattern's state bit
} pcu_coll;

void pcu_make_coll(pcu_coll* c, pcu_pattern* p, pcu_merge* m);
void pcu_begin_coll(pcu_coll* c, void* data, size_t size);
//returns false when done
bool pcu_progress_coll(pcu_coll* c);

void pcu_reduce(pcu_coll* c, pcu_merge* m, void* data, size_t size);
void pcu_bcast(pcu_coll* c, void* data, size_t size);
void pcu_allreduce(pcu_coll* c, pcu_merge* m, void* data, size_t size);
void pcu_scan(pcu_coll* c, pcu_merge* m, void* data, size_t size);

void pcu_begin_barrier(pcu_coll* c);
bool pcu_barrier_done(pcu_coll* c);
void pcu_barrier(pcu_coll* c);

#endif //PCU_COLL_H
