/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_coll.h"
#include "reel.h"
#include <string.h>

#define MIN(a,b) (((b)<(a))?(b):(a))
#define MAX(a,b) (((b)>(a))?(b):(a))

static int floor_log2(int n)
{
  int r = 0;
  while ((n >>= 1)) ++r;
  return r;
}

static int ceil_log2(int n)
{
  int r = floor_log2(n);
  if ((1 << r)<n) ++r;
  return r;
}

void pcu_merge_assign(int peers, int bit, void* local, void* incoming,
                      size_t size)
{
  (void) peers, (void) bit;
  memcpy(local,incoming,size);
}

/* initiates non-blocking calls for this
   communication step */
static void begin_coll_step(pcu_mpi_t* mpi, pcu_coll* c)
{
  int action = c->pattern->action(mpi, c->bit);
  if (action == pcu_coll_idle)
    return;
  c->message.peer = c->pattern->peer(mpi, c->bit);
  if (action == pcu_coll_send)
    pcu_mpi_send(mpi, &(c->message),mpi->coll_comm);
}

/* tries to complete this communication step.
   Returns false if communication is not done,
   otherwise wraps up communication, merges
   if necessary, and returns true */
static bool end_coll_step(pcu_mpi_t* mpi, pcu_coll* c)
{
  int action = c->pattern->action(mpi, c->bit);
  if (action == pcu_coll_idle)
    return true;
  if (action == pcu_coll_send)
    return pcu_mpi_done(mpi, &(c->message));
  pcu_message incoming;
  pcu_make_message(&incoming);
  incoming.peer = c->pattern->peer(mpi, c->bit);
  if ( ! pcu_mpi_receive(mpi, &incoming,mpi->coll_comm))
    return false;
  if (c->message.buffer.size != incoming.buffer.size)
    reel_fail("PCU unexpected incoming message.\n"
              "Most likely a PCU collective was not called by all ranks.");
  c->merge(pcu_mpi_size(mpi), c->bit, c->message.buffer.start,
           incoming.buffer.start, incoming.buffer.size);
  pcu_free_message(&incoming);
  return true;
}

void pcu_make_coll(pcu_mpi_t* mpi, pcu_coll* c, pcu_pattern* p, pcu_merge* m)
{
  // silence warning
  (void)mpi;
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
void pcu_begin_coll(pcu_mpi_t* mpi, pcu_coll* c, void* data, size_t size)
{
  pcu_set_buffer(&(c->message.buffer),data,size);
  c->bit = c->pattern->begin_bit(mpi);
  if (c->pattern->end_bit(mpi, c->bit))
    return;
  begin_coll_step(mpi, c);
}

/* makes progress on a collective operation
   started by pcu_begin_coll.
   returns false if its done. */
bool pcu_progress_coll(pcu_mpi_t* mpi, pcu_coll* c)
{
  if (c->pattern->end_bit(mpi, c->bit))
    return false;
  if (end_coll_step(mpi, c))
  {
    c->bit = c->pattern->shift(mpi, c->bit);
    if (c->pattern->end_bit(mpi, c->bit))
      return false;
    begin_coll_step(mpi, c);
  }
  return true;
}

/* reduce merges odd ranks into even ones,
   then odd multiples of 2 into even ones, etc...
   until rank 0 has all inputs merged */

static int reduce_begin_bit(pcu_mpi_t* mpi)
{
  // silence warning
  (void)mpi;
  return 1;
}

static bool reduce_end_bit(pcu_mpi_t* mpi, int bit)
{
  int rank = pcu_mpi_rank(mpi);
  if (rank==0)
    return bit >= pcu_mpi_size(mpi);
  return (bit>>1) & rank;
}

static int reduce_peer(pcu_mpi_t* mpi, int bit)
{
  return pcu_mpi_rank(mpi) ^ bit;
}

static int reduce_action(pcu_mpi_t* mpi, int bit)
{
  if (reduce_peer(mpi, bit) >= pcu_mpi_size(mpi))
    return pcu_coll_idle;
  if (bit & pcu_mpi_rank(mpi))
    return pcu_coll_send;
  return pcu_coll_recv;
}

static int reduce_shift(pcu_mpi_t* mpi, int bit)
{
  // silence warning
  (void)mpi;
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

static int bcast_begin_bit(pcu_mpi_t* mpi)
{
  int rank = pcu_mpi_rank(mpi);
  if (rank == 0)
    return 1 << ceil_log2(pcu_mpi_size(mpi));
  int bit = 1;
  while ( ! (bit & rank)) bit <<= 1;
  return bit;
}

static bool bcast_end_bit(pcu_mpi_t * mpi, int bit)
{
  // silence warning
  (void)mpi;
  return bit == 0;
}

static int bcast_peer(pcu_mpi_t * mpi, int bit)
{
  return pcu_mpi_rank(mpi) ^ bit;
}

static int bcast_action(pcu_mpi_t* mpi, int bit)
{
  if (bcast_peer(mpi, bit) >= pcu_mpi_size(mpi))
    return pcu_coll_idle;
  if (bit & pcu_mpi_rank(mpi))
    return pcu_coll_recv;
  return pcu_coll_send;
}

static int bcast_shift(pcu_mpi_t * mpi, int bit)
{
  // silence warning
  (void)mpi;
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

static int scan_up_begin_bit(pcu_mpi_t* mpi)
{
  // silence warning
  (void)mpi;
  return 1;
}

static bool scan_up_end_bit(pcu_mpi_t* mpi, int bit)
{
  return bit == (1 << floor_log2(pcu_mpi_size(mpi)));
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

static int scan_up_action(pcu_mpi_t* mpi, int bit)
{
  int rank = pcu_mpi_rank(mpi);
  if ((scan_up_could_receive(rank,bit))&&
      (0 <= scan_up_sender_for(rank,bit)))
    return pcu_coll_recv;
  int receiver = scan_up_receiver_for(rank,bit);
  if ((receiver < pcu_mpi_size(mpi))&&
      (scan_up_could_receive(receiver,bit)))
    return pcu_coll_send;
  return pcu_coll_idle;
}

static int scan_up_peer(pcu_mpi_t* mpi, int bit)
{
  int rank = pcu_mpi_rank(mpi);
  int sender = scan_up_sender_for(rank,bit);
  if ((scan_up_could_receive(rank,bit))&&
      (0 <= sender))
    return sender;
  int receiver = scan_up_receiver_for(rank,bit);
  if ((receiver < pcu_mpi_size(mpi))&&
      (scan_up_could_receive(receiver,bit)))
    return receiver;
  return -1;
}

static int scan_up_shift(pcu_mpi_t* mpi, int bit)
{
  // silence warning
  (void)mpi;
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

static int scan_down_begin_bit(pcu_mpi_t* mpi)
{
  return 1 << floor_log2(pcu_mpi_size(mpi));
}

static bool scan_down_end_bit(pcu_mpi_t* mpi, int bit)
{
  // silence warning
  (void)mpi;
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

static int scan_down_action(pcu_mpi_t * mpi, int bit)
{
  int rank = pcu_mpi_rank(mpi);
  if ((scan_down_could_send(rank,bit))&&
      (scan_down_receiver_for(rank,bit) < pcu_mpi_size(mpi)))
    return pcu_coll_send;
  int sender = scan_down_sender_for(rank,bit);
  if ((0 <= sender)&&
      (scan_down_could_send(sender,bit)))
    return pcu_coll_recv;
  return pcu_coll_idle;
}

static int scan_down_peer(pcu_mpi_t * mpi, int bit)
{
  int rank = pcu_mpi_rank(mpi);
  if (scan_down_could_send(rank,bit))
  {
    int receiver = scan_down_receiver_for(rank,bit);
    if (receiver < pcu_mpi_size(mpi))
      return receiver;
  }
  int sender = scan_down_sender_for(rank,bit);
  if ((0 <= sender)&&
      (scan_down_could_send(sender,bit)))
    return sender;
  return -1;
}

static int scan_down_shift(pcu_mpi_t* mpi, int bit)
{
  // silence warning
  (void)mpi;
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

static pcu_pattern gather =
{
  .begin_bit = reduce_begin_bit,
  .end_bit = reduce_end_bit,
  .action = reduce_action,
  .peer = reduce_peer,
  .shift = reduce_shift,
};

void pcu_reduce(pcu_mpi_t* mpi, pcu_coll* c, pcu_merge* m, void* data, size_t size)
{
  pcu_make_coll(mpi, c,&reduce,m);
  pcu_begin_coll(mpi, c,data,size);
  while(pcu_progress_coll(mpi, c));
}

void pcu_bcast(pcu_mpi_t * mpi, pcu_coll* c, void* data, size_t size)
{
  pcu_make_coll(mpi, c,&bcast,pcu_merge_assign);
  pcu_begin_coll(mpi, c,data,size);
  while(pcu_progress_coll(mpi, c));
}

void pcu_allreduce(pcu_mpi_t* mpi, pcu_coll* c, pcu_merge* m, void* data, size_t size)
{
  pcu_reduce(mpi, c,m,data,size);
  pcu_bcast(mpi, c,data,size);
}

void pcu_scan(pcu_mpi_t* mpi, pcu_coll* c, pcu_merge* m, void* data, size_t size)
{
  pcu_make_coll(mpi, c,&scan_up,m);
  pcu_begin_coll(mpi, c,data,size);
  while(pcu_progress_coll(mpi, c));
  pcu_make_coll(mpi, c,&scan_down,m);
  pcu_begin_coll(mpi, c,data,size);
  while(pcu_progress_coll(mpi, c));
}

void pcu_merge_gather(int peers, int bit, void *local, void *incoming,
                      size_t size) {
  size_t block_size = size / peers;
  // local has `bit` blocks.
  // incoming may have `bit` (if peers is a power of 2) or `bit - 1` blocks.
  // either way, writing `size - bit * block_size` prevents buffer overrun.
  // Also, all incoming blocks are from greater ranks, so they go to the right.
  memcpy(local + bit * block_size, incoming, size - bit * block_size);
}

void pcu_gather(pcu_mpi_t* mpi, pcu_coll* c, const void *send_data,
                   void *recv_data, size_t size) {
  memcpy(recv_data, send_data, size);
  pcu_make_coll(mpi, c, &gather, pcu_merge_gather);
  pcu_begin_coll(mpi, c, recv_data, size * pcu_mpi_size(mpi));
  while (pcu_progress_coll(mpi, c));
}

void pcu_allgather(pcu_mpi_t* mpi, pcu_coll* c, const void *send_data,
                   void *recv_data, size_t size) {
  pcu_gather(mpi, c, send_data, recv_data, size);
  pcu_bcast(mpi, c, recv_data, size * pcu_mpi_size(mpi));
}

/* a barrier is just an allreduce of nothing in particular */
void pcu_begin_barrier(pcu_mpi_t* mpi, pcu_coll* c)
{
  pcu_make_coll(mpi, c,&reduce,pcu_merge_assign);
  pcu_begin_coll(mpi, c,NULL,0);
}

bool pcu_barrier_done(pcu_mpi_t* mpi, pcu_coll* c)
{
  if (c->pattern == &reduce)
    if ( ! pcu_progress_coll(mpi, c))
    {
      pcu_make_coll(mpi, c,&bcast,pcu_merge_assign);
      pcu_begin_coll(mpi, c,c->message.buffer.start,c->message.buffer.size);
    }
  if (c->pattern == &bcast)
    if ( ! pcu_progress_coll(mpi,c))
      return true;
  return false;
}

void pcu_barrier(pcu_mpi_t* mpi, pcu_coll* c)
{
  // silence warning
  (void)mpi;
  pcu_begin_barrier(mpi, c);
  while( ! pcu_barrier_done(mpi, c));
}
