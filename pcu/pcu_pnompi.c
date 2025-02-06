/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_pnompi.h"
#include "pcu_buffer.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

void pcu_pmpi_send2(const pcu_mpi_t *, pcu_message* m, int tag, MPI_Comm comm);
bool pcu_pmpi_receive2(const pcu_mpi_t *, pcu_message* m, int tag, MPI_Comm);

//
// ------------------------------------------------------------------
// MPI related messages
// ------------------------------------------------------------------

//
// message as a doubled link list
typedef struct _NoMpiMsg
{
    char*             msg;      // the message
    int               size;     // message size in byte
    int               tag;      // tag associated with the message
    int               sender;   // sender rank
    int               receiver; // receiver rank
    struct _NoMpiMsg* next;     // next message
    struct _NoMpiMsg* prev;     // previous message
} NoMpiMsg;

static NoMpiMsg* _header = 0; // message headers

void add_nompi_msg(void* msg, int size, int tag, int sender, int receiver)
{
    NoMpiMsg *nmsg = malloc(sizeof(NoMpiMsg));
    char* m = malloc(sizeof(char)*size);
    char* mmsg = msg;
    for (int i=0; i<size; i++)
    {
        m[i] = mmsg[i];
    }
    nmsg->msg      = m;
    nmsg->size     = size;
    nmsg->tag      = tag;
    nmsg->sender   = sender;
    nmsg->receiver = receiver;
    nmsg->next     = 0;
    nmsg->prev     = 0;
    
    if (_header==0) {
        _header = nmsg;
        nmsg->next = nmsg;
        nmsg->prev = nmsg;
    }
    else {
        nmsg->prev = _header->prev;
        _header->prev->next = nmsg;
        _header->prev = nmsg;
        nmsg->next = _header;
    }
}

NoMpiMsg* get_nompi_msg(int tag, int receiver)
{
    NoMpiMsg* item=_header;
    if (item!=0)
    {
        do
        {
            if (item->tag == tag && item->receiver == receiver)
            {
                item->prev->next = item->next;
                item->next->prev = item->prev;
                item->next = 0;
                item->prev = 0;
                if (_header==item)
                    _header = 0;
                return item;
            }
            item = item->next;
        } while(item!=_header);
    }
    return 0;
}

void free_nompi_msg(NoMpiMsg* msg)
{
    if (msg->msg)
        free(msg->msg);
    free(msg);
}


//
// -------------------------------------------------------------------
// pcu_mpi_* functionalities
// -------------------------------------------------------------------
//

double MPI_Wtime(void)
{
  struct timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  return (double)now.tv_sec + now.tv_nsec * 1.0e-9;
}

void pcu_pmpi_init(MPI_Comm comm, pcu_mpi_t *self) {
  self->original_comm = comm;
  self->user_comm = comm+1;
  self->coll_comm = comm+2;
  self->size = 1;
  self->rank = 0;
}

void pcu_pmpi_finalize(pcu_mpi_t* self) {
  self->user_comm = 0;
  self->coll_comm = 0;
}

int pcu_pmpi_free(MPI_Comm *c) {
  (void) c;
  return 0;
}

int pcu_pmpi_split(MPI_Comm cm, int c, int k, MPI_Comm *cm2) {
  (void) cm, (void) c, (void) k, (void) cm2;
  return 1;
}

int pcu_pmpi_size(const pcu_mpi_t *self) {
  return self->size;
}

int pcu_pmpi_rank(const pcu_mpi_t *self) {
  return self->rank;
}

void pcu_pmpi_send(const pcu_mpi_t *self, pcu_message* m, MPI_Comm comm) {
  pcu_pmpi_send2(self, m, 0, comm);
}

void pcu_pmpi_send2(const pcu_mpi_t *self, pcu_message* m, int tag,
  MPI_Comm c) {
  (void) c;
  if (m->buffer.size > (size_t)INT_MAX) {
    fprintf(stderr, "ERROR PCU message size exceeds INT_MAX... exiting\n");
    abort();
  }

  add_nompi_msg(
    m->buffer.start, (int)(m->buffer.size), tag, self->rank, m->peer
  );
}

bool pcu_pmpi_done(const pcu_mpi_t* a, pcu_message* b) {
  (void) a, (void) b;
  return true;
}

bool pcu_pmpi_receive(const pcu_mpi_t *self, pcu_message* m, MPI_Comm comm) {
  return pcu_pmpi_receive2(self, m, 0, comm);
}

bool pcu_pmpi_receive2(const pcu_mpi_t *self, pcu_message* m, int tag,
  MPI_Comm comm) {
  (void) comm;
  NoMpiMsg* msg = get_nompi_msg(tag, self->rank);
  if (msg==0) return false;
  m->peer = msg->sender;
  int msize = msg->size;
  pcu_resize_buffer(&(m->buffer),(size_t)msize);
  for (int i = 0; i < msize; i++) {
    ((char*)m->buffer.start)[i] = msg->msg[i];
  }
  // delete the message
  free(msg->msg);
  free(msg);
  return true;
}

int pcu_pmpi_allreduce(const void* a, void* b, int c,
  MPI_Datatype d, MPI_Op e, MPI_Comm f) {
  (void) a, (void) b, (void) c, (void) d, (void) e, (void) f;
  return 0;
}

int  pcu_pmpi_allgather(const void * a, int b, MPI_Datatype c, void * d, int e,
  MPI_Datatype f, MPI_Comm g) {
  (void) a, (void) b, (void) c, (void) d, (void) e, (void) f, (void) g;
  return 0;
}

int pcu_pmpi_barrier(MPI_Comm a) {
  (void) a;
  return 0;
}
