/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_tmpi.h"
#include "pcu_thread.h"
#include "pcu_pmpi.h"
#include "pcu_common.h"
#include <assert.h>
#include <limits.h>

pcu_mpi pcu_tmpi =
{ .size = pcu_tmpi_size,
  .rank = pcu_tmpi_rank,
  .send = pcu_tmpi_send,
  .done = pcu_tmpi_done,
  .receive = pcu_tmpi_receive };

int pcu_tmpi_size()
{
  return pcu_pmpi_size() * pcu_thread_size();
}

int pcu_tmpi_rank()
{
  return pcu_thread_size() * pcu_pmpi_rank() + pcu_thread_rank();
}

#define THREAD_BITS 10
#define TO_SHIFT 0
#define TO_WIDTH THREAD_BITS
#define FROM_SHIFT (TO_SHIFT+TO_WIDTH)
#define FROM_WIDTH THREAD_BITS

static int width_mask(int width)
{
  return (1 << width) - 1;
}

static int compose(int value, int width, int shift)
{
  return (value & width_mask(width)) << shift;
}

static int decompose(int value, int width, int shift)
{
  return (value >> shift) & width_mask(width);
}

static int from_thread_of(int compound_tag)
{
  return decompose(compound_tag,FROM_WIDTH,FROM_SHIFT);
}

static int to_thread_of(int compound_tag)
{
  return decompose(compound_tag,TO_WIDTH,TO_SHIFT);
}

static int make_compound_tag(int from_thread, int to_thread)
{
  int result = 0;
  result |= compose(to_thread,TO_WIDTH,TO_SHIFT);
  result |= compose(from_thread,FROM_WIDTH,FROM_SHIFT);
  assert(from_thread_of(result) == from_thread);
  assert(to_thread_of(result) == to_thread);
  return result;
}

void pcu_tmpi_send(pcu_message* m, MPI_Comm comm)
{
  int thread_size = pcu_thread_size();
  int thread_rank = pcu_thread_rank();
  int peer_thread = m->peer % thread_size;
  int peer_process = m->peer / thread_size;
  int tag = make_compound_tag(thread_rank,peer_thread);
  m->peer = peer_process;
  pcu_pmpi_send2(m,tag,comm);
  m->peer = peer_process*thread_size + peer_thread;
}

bool pcu_tmpi_done(pcu_message* m)
{
  return pcu_pmpi_done(m);
}

bool pcu_tmpi_receive(pcu_message* m, MPI_Comm comm)
{
  MPI_Status status;
  int flag;
  int thread_size = pcu_thread_size();
  int thread_rank = pcu_thread_rank();
  int peer_process = INT_MIN;
  int peer_thread = INT_MIN;
  int mpi_tag;
  if (m->peer == MPI_ANY_SOURCE)
    mpi_tag = MPI_ANY_TAG;
  else
  {
    peer_process = m->peer / thread_size;
    peer_thread = m->peer % thread_size;
    m->peer = peer_process;
    mpi_tag = make_compound_tag(peer_thread,thread_rank);
  }
  MPI_Iprobe(m->peer,mpi_tag,comm,&flag,&status);
  if (!flag)
  {
    if (m->peer != MPI_ANY_SOURCE)
      m->peer = peer_process * thread_size + peer_thread;
    return false;
  }
  int count;
  MPI_Get_count(&status,MPI_BYTE,&count);
  if (m->peer == MPI_ANY_SOURCE)
  {
    if ((to_thread_of(status.MPI_TAG) != thread_rank))
      return false;
    mpi_tag = status.MPI_TAG;
    peer_process = status.MPI_SOURCE;
    peer_thread = from_thread_of(mpi_tag);
    m->peer = peer_process;
  }
  pcu_resize_buffer(&(m->buffer),(size_t)count);
  MPI_Recv(
      m->buffer.start,
      (int)m->buffer.size,
      MPI_BYTE,
      m->peer,
      mpi_tag,
      comm,
      MPI_STATUS_IGNORE);
  m->peer = peer_process*thread_size + peer_thread;
  return true;
}

void pcu_tmpi_check_support(void)
{
  int provided;
  MPI_Query_thread(&provided);
  if (provided != MPI_THREAD_MULTIPLE)
    pcu_fail("MPI_Init_thread was not called with MPI_THREAD_MULTIPLE");
}
