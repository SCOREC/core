/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_msg.h"
#include "pcu_common.h"
#include "pcu_pmpi.h"
#include <string.h>

/* the pcu_msg algorithm for a communication phase
   is as follows:

1  barrier
2  pack data to be sent
3  requests = send all packed data
4  while (requests not done)
5    receive and process data
6  begin barrier
7  while (barrier not done)
8    receive and process data

   The main goal of this is to detect when no more
   messages will be received.

   To prove that this works, note the following:
   1. messages are sent using a synchronous
      non-blocking send: requests are not done
      until the message is received at its destination.
   2. the barrier is not done on any rank until
      after all ranks have begun the barrier
   3. all requests of this rank are done before this
      rank begins the barrier
   From this we can prove:
   1. if all requests of this rank are done, all messages sent by
      this rank have been received at their destinations
   2. when the barrier is done, all messages sent by all ranks
      have been received at their destinations, i.e.
      no more messages can be received.

   Finally, the barrier at line 1 is there because it takes
   different amounts of time for each rank to be notified
   that the barrier is done, and in that time the rank is
   naively processing everything it receives.
   If another rank is notified first and quickly goes on to
   a new phase, it may be able to send a message that is
   received by the slow rank out-of-phase.
*/

//enumeration for pcu_msg.state
enum {
  idle_state, //in between phases
  pack_state, //after phase start, before sending
  send_recv_state, //starting to receive, sends still going
  recv_state //sends are done, still receiving
};

static void make_comm(pcu_msg* m)
{
  pcu_make_aa(&(m->peers));
  pcu_make_message(&(m->received));
  m->state = idle_state;
}

void pcu_make_msg(pcu_msg* m)
{
  make_comm(m);
  m->file = NULL;
}

static void free_peers(pcu_aa_tree* t)
{
  if (pcu_aa_empty(*t))
    return;
  free_peers(&((*t)->left));
  free_peers(&((*t)->right));
  pcu_msg_peer* peer;
  peer = (pcu_msg_peer*) *t;
  pcu_free_message(&(peer->message));
  pcu_free(peer);
  pcu_make_aa(t);
}

void pcu_msg_start(pcu_msg* m)
{
  if (m->state != idle_state) pcu_fail("Start called at the wrong time");
  /* this barrier ensures no one starts a new superstep
     while others are receiving in the past superstep.
     It is the only blocking call in the pcu_msg system. */
  pcu_barrier(&(m->coll));
  m->state = pack_state;
}

static bool peer_less(pcu_aa_node* a, pcu_aa_node* b)
{
  return ((pcu_msg_peer*)a)->message.peer
       < ((pcu_msg_peer*)b)->message.peer;
}

static pcu_msg_peer* find_peer(pcu_aa_tree t, int id)
{
  pcu_msg_peer key;
  key.message.peer = id;
  return (pcu_msg_peer*) pcu_aa_find(&(key.node),t,peer_less);
}

static pcu_msg_peer* make_peer(int id)
{
  pcu_msg_peer* p;
  PCU_MALLOC(p,1);
  pcu_make_message(&(p->message));
  p->message.peer = id;
  return p;
}

void* pcu_msg_pack(pcu_msg* m, int id, size_t size)
{
  if (m->state != pack_state)
    pcu_fail("Pack or Write called at the wrong time");
  pcu_msg_peer* peer = find_peer(m->peers,id);
  if (!peer)
  {
    peer = make_peer(id);
    pcu_aa_insert(&(peer->node),&(m->peers),peer_less);
  }
  return pcu_push_buffer(&(peer->message.buffer),size);
}

size_t pcu_msg_packed(pcu_msg* m, int id)
{
  if (m->state != pack_state) pcu_fail("Packed called at the wrong time");
  pcu_msg_peer* peer = find_peer(m->peers,id);
  if (!peer) pcu_fail("pcu_msg_packed called but nothing was packed");
  return peer->message.buffer.size;
}

static void send_peers(pcu_aa_tree t)
{
  if (pcu_aa_empty(t))
    return;
  pcu_msg_peer* peer;
  peer = (pcu_msg_peer*)t;
  pcu_mpi_send(&(peer->message),pcu_user_comm);
  send_peers(t->left);
  send_peers(t->right);
}

void pcu_msg_send(pcu_msg* m)
{
  if (m->state != pack_state)
    pcu_fail("Send called at the wrong time");
  send_peers(m->peers);
  m->state = send_recv_state;
}

static bool done_sending_peers(pcu_aa_tree t)
{
  if (pcu_aa_empty(t))
    return true;
  pcu_msg_peer* peer;
  peer = (pcu_msg_peer*)t;
  return pcu_mpi_done(&(peer->message))
    && done_sending_peers(t->left)
    && done_sending_peers(t->right);
}

static bool receive_global(pcu_msg* m)
{
  m->received.peer = MPI_ANY_SOURCE;
  while ( ! pcu_mpi_receive(&(m->received),pcu_user_comm))
  {
    if (m->state == send_recv_state)
      if (done_sending_peers(m->peers))
      {
        pcu_begin_barrier(&(m->coll));
        m->state = recv_state;
      }
    if (m->state == recv_state)
      if (pcu_barrier_done(&(m->coll)))
        return false;
  }
  return true;
}

static void free_comm(pcu_msg* m)
{
  free_peers(&(m->peers));
  pcu_free_message(&(m->received));
}

bool pcu_msg_receive(pcu_msg* m)
{
  if ((m->state != send_recv_state)&&
      (m->state != recv_state))
    pcu_fail("Receive called at the wrong time");
  if ( ! pcu_msg_unpacked(m))
    pcu_fail("Receive called before previous message unpacked");
  if (receive_global(m))
  {
    pcu_begin_buffer(&(m->received.buffer));
    return true;
  }
  m->state = idle_state;
  free_comm(m);
  make_comm(m);
  return false;
}

void* pcu_msg_unpack(pcu_msg* m, size_t size)
{
  return pcu_walk_buffer(&(m->received.buffer),size);
}

bool pcu_msg_unpacked(pcu_msg* m)
{
  return pcu_buffer_walked(&(m->received.buffer));
}

int pcu_msg_received_from(pcu_msg* m)
{
  return m->received.peer;
}

size_t pcu_msg_received_size(pcu_msg* m)
{
  return m->received.buffer.capacity;
}

void pcu_free_msg(pcu_msg* m)
{
  free_comm(m);
  if (m->file)
    fclose(m->file);
}

