/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_coll.h"
#include "pcu_pmpi.h"
#include "pcu_common.h"
#include <assert.h>
#include <string.h>

void pcu_merge_assign(void* local, void* incoming, size_t size)
{
  memcpy(local,incoming,size);
}

void pcu_add_doubles(void* local, void* incoming, size_t size)
{
  double* a = local;
  double* b= incoming;
  size_t n = size/sizeof(double);
  for (size_t i=0; i < n; ++i)
    a[i] += b[i];
}

void pcu_max_doubles(void* local, void* incoming, size_t size)
{
  double* a = local;
  double* b= incoming;
  size_t n = size/sizeof(double);
  for (size_t i=0; i < n; ++i)
    a[i] = MAX(a[i],b[i]);
}

void pcu_min_doubles(void* local, void* incoming, size_t size)
{
  double* a = local;
  double* b= incoming;
  size_t n = size/sizeof(double);
  for (size_t i=0; i < n; ++i)
    a[i] = MIN(a[i],b[i]);
}

void pcu_add_ints(void* local, void* incoming, size_t size)
{
  int* a = local;
  int* b= incoming;
  size_t n = size/sizeof(int);
  for (size_t i=0; i < n; ++i)
    a[i] += b[i];
}

void pcu_min_ints(void* local, void* incoming, size_t size)
{
  int* a = local;
  int* b= incoming;
  size_t n = size/sizeof(int);
  for (size_t i=0; i < n; ++i)
    a[i] = MIN(a[i],b[i]);
}

void pcu_max_ints(void* local, void* incoming, size_t size)
{
  int* a = local;
  int* b= incoming;
  size_t n = size/sizeof(int);
  for (size_t i=0; i < n; ++i)
    a[i] = MAX(a[i],b[i]);
}

void pcu_add_longs(void* local, void* incoming, size_t size)
{
  long* a = local;
  long* b= incoming;
  size_t n = size/sizeof(long);
  for (size_t i=0; i < n; ++i)
    a[i] += b[i];
}

/* initiates non-blocking calls for this
   communication step */
static void begin_coll_step(pcu_coll* c)
{
  int action = c->pattern->action(c->bit);
  if (action == pcu_coll_idle)
    return;
  c->message.peer = c->pattern->peer(c->bit);
  if (action == pcu_coll_send)
    pcu_mpi_send(&(c->message),pcu_coll_comm);
}

/* tries to complete this communication step.
   Returns false if communication is not done,
   otherwise wraps up communication, merges
   if necessary, and returns true */
static bool end_coll_step(pcu_coll* c)
{
  int action = c->pattern->action(c->bit);
  if (action == pcu_coll_idle)
    return true;
  if (action == pcu_coll_send)
    return pcu_mpi_done(&(c->message));
  pcu_message incoming;
  pcu_make_message(&incoming);
  incoming.peer = c->pattern->peer(c->bit);
  if ( ! pcu_mpi_receive(&incoming,pcu_coll_comm))
    return false;
  if (c->message.buffer.size != incoming.buffer.size)
    pcu_fail("collective not called by all ranks or pcu bug");
  c->merge(c->message.buffer.start,incoming.buffer.start,incoming.buffer.size);
  pcu_free_message(&incoming);
  return true;
}

void pcu_make_coll(pcu_coll* c, pcu_pattern* p, pcu_merge* m)
{
  c->pattern = p;
  c->merge = m;
}

/* the abstract algorithm for a collective communication
   pattern is as follows:
   for (bit = begin_bit; ! end_bit(bit); bit = shift(bit))
     begin_step(bit)
     end_step(bit)

   to run this as a non-blocking operation, we split
   it between the begin/end steps.
   pcu_coll will return "not yet done" after begin_step
   and will pick up at end_step when called to progress again

   this becomes messy considering in some cases begin_step
   is not called at all, see pcu_begin_coll and pcu_progress_coll
   below:
*/

/* begins a non-blocking collective.
   The collective operation should be set with pcu_make_coll first.
   data[0..size] is the input/output local data */
void pcu_begin_coll(pcu_coll* c, void* data, size_t size)
{
  pcu_set_buffer(&(c->message.buffer),data,size);
  c->bit = c->pattern->begin_bit();
  if (c->pattern->end_bit(c->bit))
    return;
  begin_coll_step(c);
}

/* makes progress on a collective operation
   started by pcu_begin_coll.
   returns false if its done. */
bool pcu_progress_coll(pcu_coll* c)
{
  if (c->pattern->end_bit(c->bit))
    return false;
  if (end_coll_step(c))
  {
    c->bit = c->pattern->shift(c->bit);
    if (c->pattern->end_bit(c->bit))
      return false;
    begin_coll_step(c);
  }
  return true;
}

/* reduce merges odd ranks into even ones,
   then odd multiples of 2 into even ones, etc...
   until rank 0 has all inputs merged */

static int reduce_begin_bit(void)
{
  return 1;
}

static bool reduce_end_bit(int bit)
{
  int rank = pcu_mpi_rank();
  if (rank==0)
    return bit >= pcu_mpi_size();
  return (bit>>1) & rank;
}

static int reduce_peer(int bit)
{
  return pcu_mpi_rank() ^ bit;
}

static int reduce_action(int bit)
{
  if (reduce_peer(bit) >= pcu_mpi_size())
    return pcu_coll_idle;
  if (bit & pcu_mpi_rank())
    return pcu_coll_send;
  return pcu_coll_recv;
}

static int reduce_shift(int bit)
{
  return bit << 1;
}

static pcu_pattern reduce =
{
  .begin_bit = reduce_begin_bit,
  .end_bit = reduce_end_bit,
  .action = reduce_action,
  .peer = reduce_peer,
  .shift = reduce_shift,
};

/* broadcast is the opposite of reduce,
   the pattern runs backwards and send/recv
   are flipped. */

static int bcast_begin_bit(void)
{
  int rank = pcu_mpi_rank();
  if (rank == 0)
    return 1 << pcu_ceil_log2(pcu_mpi_size());
  int bit = 1;
  while ( ! (bit & rank)) bit <<= 1;
  return bit;
}

static bool bcast_end_bit(int bit)
{
  return bit == 0;
}

static int bcast_peer(int bit)
{
  return pcu_mpi_rank() ^ bit;
}

static int bcast_action(int bit)
{
  if (bcast_peer(bit) >= pcu_mpi_size())
    return pcu_coll_idle;
  if (bit & pcu_mpi_rank())
    return pcu_coll_recv;
  return pcu_coll_send;
}

static int bcast_shift(int bit)
{
  return bit >> 1;
}

static pcu_pattern bcast =
{
  .begin_bit = bcast_begin_bit,
  .end_bit = bcast_end_bit,
  .action = bcast_action,
  .peer = bcast_peer,
  .shift = bcast_shift,
};

/* The scan patterns are from the first
   algorithm reviewed by this paper:

   Peter Sanders, Jesper Larsson Traff.
   "Parallel Prefix (Scan) Algorithms for MPI".
*/

static int scan_up_begin_bit(void)
{
  return 1;
}

static bool scan_up_end_bit(int bit)
{
  return bit == (1 << pcu_floor_log2(pcu_mpi_size()));
}

static bool scan_up_could_receive(int rank, int bit)
{
  int mask = ((bit << 1)-1);
  return ((rank & mask) == mask);
}

static int scan_up_sender_for(int rank, int bit)
{
  return rank - bit;
}

static int scan_up_receiver_for(int rank, int bit)
{
  return rank + bit;
}

static int scan_up_action(int bit)
{
  int rank = pcu_mpi_rank();
  if ((scan_up_could_receive(rank,bit))&&
      (0 <= scan_up_sender_for(rank,bit)))
    return pcu_coll_recv;
  int receiver = scan_up_receiver_for(rank,bit);
  if ((receiver < pcu_mpi_size())&&
      (scan_up_could_receive(receiver,bit)))
    return pcu_coll_send;
  return pcu_coll_idle;
}

static int scan_up_peer(int bit)
{
  int rank = pcu_mpi_rank();
  int sender = scan_up_sender_for(rank,bit);
  if ((scan_up_could_receive(rank,bit))&&
      (0 <= sender))
    return sender;
  int receiver = scan_up_receiver_for(rank,bit);
  if ((receiver < pcu_mpi_size())&&
      (scan_up_could_receive(receiver,bit)))
    return receiver;
  return -1;
}

static int scan_up_shift(int bit)
{
  return bit << 1;
}

static pcu_pattern scan_up =
{
  .begin_bit = scan_up_begin_bit,
  .end_bit = scan_up_end_bit,
  .action = scan_up_action,
  .peer = scan_up_peer,
  .shift = scan_up_shift,
};

static int scan_down_begin_bit(void)
{
  return 1 << pcu_floor_log2(pcu_mpi_size());
}

static bool scan_down_end_bit(int bit)
{
  return bit == 1;
}

static bool scan_down_could_send(int rank, int bit)
{
  int mask = (bit-1);
  return ((rank & mask) == mask);
}

static int scan_down_receiver_for(int rank, int bit)
{
  return rank + (bit >> 1);
}

static int scan_down_sender_for(int rank, int bit)
{
  return rank - (bit >> 1);
}

static int scan_down_action(int bit)
{
  int rank = pcu_mpi_rank();
  if ((scan_down_could_send(rank,bit))&&
      (scan_down_receiver_for(rank,bit) < pcu_mpi_size()))
    return pcu_coll_send;
  int sender = scan_down_sender_for(rank,bit);
  if ((0 <= sender)&&
      (scan_down_could_send(sender,bit)))
    return pcu_coll_recv;
  return pcu_coll_idle;
}

static int scan_down_peer(int bit)
{
  int rank = pcu_mpi_rank();
  if (scan_down_could_send(rank,bit))
  {
    int receiver = scan_down_receiver_for(rank,bit);
    if (receiver < pcu_mpi_size())
      return receiver;
  }
  int sender = scan_down_sender_for(rank,bit);
  if ((0 <= sender)&&
      (scan_down_could_send(sender,bit)))
    return sender;
  return -1;
}

static int scan_down_shift(int bit)
{
  return bit >> 1;
}

static pcu_pattern scan_down =
{
  .begin_bit = scan_down_begin_bit,
  .end_bit = scan_down_end_bit,
  .action = scan_down_action,
  .peer = scan_down_peer,
  .shift = scan_down_shift,
};

void pcu_reduce(pcu_coll* c, pcu_merge* m, void* data, size_t size)
{
  pcu_make_coll(c,&reduce,m);
  pcu_begin_coll(c,data,size);
  while(pcu_progress_coll(c));
}

void pcu_bcast(pcu_coll* c, void* data, size_t size)
{
  pcu_make_coll(c,&bcast,pcu_merge_assign);
  pcu_begin_coll(c,data,size);
  while(pcu_progress_coll(c));
}

void pcu_allreduce(pcu_coll* c, pcu_merge* m, void* data, size_t size)
{
  pcu_reduce(c,m,data,size);
  pcu_bcast(c,data,size);
}

void pcu_scan(pcu_coll* c, pcu_merge* m, void* data, size_t size)
{
  pcu_make_coll(c,&scan_up,m);
  pcu_begin_coll(c,data,size);
  while(pcu_progress_coll(c));
  pcu_make_coll(c,&scan_down,m);
  pcu_begin_coll(c,data,size);
  while(pcu_progress_coll(c));
}

/* a barrier is just an allreduce of nothing in particular */
void pcu_begin_barrier(pcu_coll* c)
{
  pcu_make_coll(c,&reduce,pcu_merge_assign);
  pcu_begin_coll(c,NULL,0);
}

bool pcu_barrier_done(pcu_coll* c)
{
  if (c->pattern == &reduce)
    if ( ! pcu_progress_coll(c))
    {
      pcu_make_coll(c,&bcast,pcu_merge_assign);
      pcu_begin_coll(c,c->message.buffer.start,c->message.buffer.size);
    }
  if (c->pattern == &bcast)
    if ( ! pcu_progress_coll(c))
      return true;
  return false;
}

void pcu_barrier(pcu_coll* c)
{
  pcu_begin_barrier(c);
  while( ! pcu_barrier_done(c));
}
