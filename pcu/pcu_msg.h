/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_MSG_H
#define PCU_MSG_H

#include "pcu_coll.h"
#include "pcu_aa.h"
#include "pcu_io.h"

/* the PCU Messenger (pcu_msg for short) system implements
   a non-blocking Bulk Synchronous Parallel communication model
   it is based on PCU non-blocking Collectives and pcu_mpi,
   so it also works in hybrid mode */

/* a structure containing all data to be sent to a peer during
   this communication phase. */
typedef struct
{
  pcu_aa_node node; //binary tree node for lookup
  pcu_message message; //send buffer and peer id
} pcu_msg_peer;

struct pcu_order_struct;

struct pcu_msg_struct
{
  pcu_aa_tree peers; //binary tree of send buffers
  pcu_message received; //current received buffer
  pcu_coll coll; //collective operation object
  int state; //state within a communication phase
  /* below this point are variables that just need
     to be thread-specific but have been tacked onto
     pcu_msg. if this gets out of hand, create a
     pcu_thread struct to or something */
  FILE* file; //messenger-unique input or output file
  struct pcu_order_struct* order;
};
typedef struct pcu_msg_struct pcu_msg;

void pcu_make_msg(pcu_msg* m);
void pcu_msg_start(pcu_msg* b);
void* pcu_msg_pack(pcu_msg* m, int id, size_t size);
#define PCU_MSG_PACK(m,id,o) \
memcpy(pcu_msg_pack(m,id,sizeof(o)),&(o),sizeof(o))
size_t pcu_msg_packed(pcu_msg* m, int id);
void pcu_msg_send(pcu_msg* m);
bool pcu_msg_receive(pcu_msg* m);
void* pcu_msg_unpack(pcu_msg* m, size_t size);
#define PCU_MSG_UNPACK(m,o) \
memcpy(&(o),pcu_msg_unpack(m,sizeof(o)),sizeof(o))
bool pcu_msg_unpacked(pcu_msg* m);
int pcu_msg_received_from(pcu_msg* m);
size_t pcu_msg_received_size(pcu_msg* m);
void pcu_free_msg(pcu_msg* m);

#endif //PCU_MSG_H
