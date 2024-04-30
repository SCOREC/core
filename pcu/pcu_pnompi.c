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

static int global_size;
static int global_rank;

MPI_Comm original_comm;
MPI_Comm pcu_user_comm;
MPI_Comm pcu_coll_comm;


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
    nmsg->receiver = sender;
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
  return 0.0;
}

pcu_mpi pcu_pmpi =
{ .size = pcu_pmpi_size,
  .rank = pcu_pmpi_rank,
  .send = pcu_pmpi_send,
  .done = pcu_pmpi_done,
  .receive = pcu_pmpi_receive };

void pcu_pmpi_init(MPI_Comm comm)
{
    original_comm = comm;
    pcu_user_comm = comm+1;
    pcu_coll_comm = comm+2;

    global_size = 1;
    global_rank = 0;
    
//  MPI_Comm_dup(comm,&pcu_user_comm);
//  MPI_Comm_dup(comm,&pcu_coll_comm);
//  MPI_Comm_size(comm,&global_size);
//  MPI_Comm_rank(comm,&global_rank);
}

void pcu_pmpi_finalize(void)
{
    pcu_user_comm = 0;
    pcu_coll_comm = 0;
//  MPI_Comm_free(&pcu_user_comm);
//  MPI_Comm_free(&pcu_coll_comm);
}

int pcu_pmpi_free(MPI_Comm* comm)
{
  return 0;
}

int pcu_pmpi_split(MPI_Comm c, int color, int key, MPI_Comm* nc)
{
  return 1;
}

int pcu_pmpi_size(void)
{
  return global_size;
}

int pcu_pmpi_rank(void)
{
  return global_rank;
}

void pcu_pmpi_send(pcu_message* m, MPI_Comm comm)
{
  pcu_pmpi_send2(m,0,comm);
}

void pcu_pmpi_send2(pcu_message* m, int tag, MPI_Comm comm)
{
    if( m->buffer.size > (size_t)INT_MAX ) {
        fprintf(stderr, "ERROR PCU message size exceeds INT_MAX... exiting\n");
        abort();
    }
    
    add_nompi_msg(m->buffer.start,
                  (int)(m->buffer.size),
                  tag,
                  global_rank,
                  m->peer
                  );
}

bool pcu_pmpi_done(pcu_message* m)
{
    return true;
//  int flag;
//  MPI_Test(&(m->request),&flag,MPI_STATUS_IGNORE);
//  return flag;
}

bool pcu_pmpi_receive(pcu_message* m, MPI_Comm comm)
{
    return pcu_pmpi_receive2(m,0,comm);
}

bool pcu_pmpi_receive2(pcu_message* m, int tag, MPI_Comm comm)
{
    NoMpiMsg* msg = get_nompi_msg(tag,global_rank);
    if (msg==0)
    {
        printf("Did not receive a proper message\n");
        return false;
    }
    
    int msize = msg->size;
    
    pcu_resize_buffer(&(m->buffer),(size_t)msize);
    
    for (int i=0; i<msize; i++)
    {
        ((char*)m->buffer.start)[i] = msg->msg[i];
    }
    
    // delete the message
    free(msg->msg);
    free(msg);
    return true;
}

void pcu_pmpi_switch(MPI_Comm new_comm)
{
    return ;
//  pcu_pmpi_finalize();
//  pcu_pmpi_init(new_comm);
}

MPI_Comm pcu_pmpi_comm(void)
{
  return original_comm;
}

int pcu_pmpi_allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return 0;
}

int  pcu_pmpi_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  return 0;
}

int  pcu_pmpi_barrier(MPI_Comm comm)
{
  return 0;
}
