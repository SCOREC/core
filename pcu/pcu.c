/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
/** \file pcu.c
    \brief The PCU communication interface */
/** \page pcu PCU
  PCU (the Parallel Control Utility) is a library for parallel computation
  based on MPI with additional support for hybrid MPI/thread environments.
  PCU provides three things to users:
    1. A hybrid phased message passing system
    2. Hybrid collective operations
    3. A thread management system

  Phased message passing is similar to Bulk Synchronous Parallel.
  All messages are exchanged in a phase, which is a collective operation
  involving all threads in the parallel program.
  During a phase, the following events happen in sequence:
    1. All threads send non-blocking messages to other threads
    2. All threads receive all messages sent to them during this phase
  PCU provides termination detection, which is the ability to detect when all
  messages have been received without prior knowledge of which threads
  are sending to which.

  To write hybrid MPI/thread programs, PCU provides a function that creates
  threads within an MPI process, similar to the way mpirun creates multiple
  processes. PCU assigns ranks to these threads and has them each run the same
  function, with thread-specific input arguments to the function.

  Once a program has created threads using PCU, it can call the message passing
  API from within threads, which will behave as if each thread
  were an MPI process. Threads have unique ranks and can send messages
  to one another, regardless of which process they are in.

  The API documentation is here: pcu.c
*/

#include <string.h>
#include <stdarg.h>
#include "PCUConfig.h"
#include "PCU.h"
#include "pcu_common.h"
#include "pcu_msg.h"
#include "pcu_pmpi.h"
#include "pcu_order.h"

#if ENABLE_THREADS
#include "pcu_thread.h"
#include "pcu_tmpi.h"
#else
static void fail_no_threads(void) __attribute__((noreturn));
static void fail_no_threads(void)
{
  pcu_fail("threads disabled. reconfigure with --enable-threads");
}
#endif

enum state { uninit, init };
static enum state global_state = uninit;
static pcu_msg global_pmsg;
#if ENABLE_THREADS
static pcu_msg* global_tmsg = NULL;
static PCU_Thrd_Func global_function = NULL;
static void** global_args = NULL;
#endif

static pcu_msg* get_msg()
{
#if ENABLE_THREADS
  if (pcu_get_mpi() == &pcu_tmpi)
    return global_tmsg + pcu_thread_rank();
#endif
  return &global_pmsg;
}

/** \brief Initializes the PCU library.
  \details This function must be called by all MPI processes before
  calling any other PCU functions.
  MPI_Init or MPI_Init_thread should be called before this function.
 */
int PCU_Comm_Init(void)
{
  if (global_state != uninit)
    pcu_fail("nested calls to Comm_Init");
  pcu_pmpi_init(MPI_COMM_WORLD);
  pcu_set_mpi(&pcu_pmpi);
  pcu_make_msg(&global_pmsg);
  global_state = init;
  return PCU_SUCCESS;
}

/** \brief Frees all PCU library structures.
  \details This function must be called by all MPI processes after all
  other calls to PCU, and before calling MPI_Finalize.
 */
int PCU_Comm_Free(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Free called before Comm_Init");
  if (global_pmsg.order)
    pcu_order_free(global_pmsg.order);
  pcu_free_msg(&global_pmsg);
  pcu_pmpi_finalize();
  global_state = uninit;
  return PCU_SUCCESS;
}

/** \brief Returns the communication rank of the calling thread.
  \details when called from a non-threaded MPI process, this function is
  equivalent to MPI_Comm_rank(MPI_COMM_WORLD,rank).

  When called from a thread inside PCU_Thrd_Run, the rank is unique to a thread
  in the whole MPI job.
  Ranks are consecutive from 0 to \f$pt-1\f$ for a program with
  \f$p\f$ processes and \f$t\f$ threads per process.
  Ranks are contiguous within a process, so that the \f$t\f$ threads in process
  \f$i\f$ are numbered from \f$ti\f$ to \f$ti+t-1\f$.
 */
int PCU_Comm_Self(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Self called before Comm_Init");
  return pcu_mpi_rank();
}

/** \brief Returns the number of threads in the program.
  \details when called from a non-threaded MPI process, this function is
  equivalent to MPI_Comm_size(MPI_COMM_WORLD,size).

  When called from a thread inside PCU_Thrd_Run, the size is \f$pt\f$, where
  \f$p\f$ is the number of MPI processes and \f$t\f$ is the number of threads
  per process, which is the nthreads argument passed to PCU_Thrd_Run.
 */
int PCU_Comm_Peers(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Peers called before Comm_Init");
  return pcu_mpi_size();
}

/** \brief Begins a PCU communication phase.
  \details This function must be called by all threads in the MPI job
  at the beginning of each phase of communication.
  After calling this function, each thread may call functions like
  PCU_Comm_Pack or PCU_Comm_Write.
*/
void PCU_Comm_Begin(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Begin called before Comm_Init");
  pcu_msg_start(get_msg());
}

/** \brief Packs data to be sent to \a to_rank.
  \details This function appends the block of \a size bytes starting
  at \a data to the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Pack(int to_rank, const void* data, size_t size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Pack called before Comm_Init");
  if ((to_rank < 0)||(to_rank >= pcu_mpi_size()))
    pcu_fail("Invalid rank in Comm_Pack");
  memcpy(pcu_msg_pack(get_msg(),to_rank,size),data,size);
  return PCU_SUCCESS;
}

/** \brief Sends all buffers for this communication phase.
  \details This function should be called by all threads in the MPI job
  after calls to PCU_Comm_Pack or PCU_Comm_Write and before calls
  to PCU_Comm_Listen or PCU_Comm_Read.
  All buffers from this thread are sent out and receiving
  may begin after this call.
 */
int PCU_Comm_Send(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Send called before Comm_Init");
  pcu_msg_send(get_msg());
  return PCU_SUCCESS;
}

/** \brief Tries to receive a buffer for this communication phase.
  \details Either this function or PCU_Comm_Read should be called at least
  once by all threads during the communication phase, after PCU_Comm_Send
  is called. The result will be false if and only if the communication phase
  is over and there are no more buffers to receive.
  Otherwise, a buffer was received.
  Its contents are retrievable through PCU_Comm_Unpack, and its metadata through
  PCU_Comm_Sender and PCU_Comm_Received.
  Users should unpack all data from this buffer before calling this function
  again, because the previously received buffer is destroyed by the call.
 */
bool PCU_Comm_Listen(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Listen called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    return pcu_order_receive(m->order, m);
  return pcu_msg_receive(m);
}

/** \brief Returns in * \a from_rank the sender of the current received buffer.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
int PCU_Comm_Sender(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Sender called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    return pcu_order_received_from(m->order);
  return pcu_msg_received_from(m);
}

/** \brief Returns true if the current received buffer has been unpacked.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
bool PCU_Comm_Unpacked(void)
{
  if (global_state == uninit)
    pcu_fail("Comm_Unpacked called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    return pcu_order_unpacked(m->order);
  return pcu_msg_unpacked(m);
}

/** \brief Unpacks a block of data from the current received buffer.
  \details This function should be called after a successful PCU_Comm_Listen.
  \a data must point to a block of memory of at least \a size bytes, into
  which the next \a size bytes of the current received buffer will be written.
  Subsequent calls will begin unpacking where this call left off,
  so that the entire received buffer can be unpacked by a sequence of calls to
  this function.
  Users must ensure that there remain \a size bytes to be unpacked,
  PCU_Comm_Unpacked can help with this.
 */
int PCU_Comm_Unpack(void* data, size_t size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Unpack called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    memcpy(data,pcu_order_unpack(m->order,size),size);
  else
    memcpy(data,pcu_msg_unpack(m,size),size);
  return PCU_SUCCESS;
}

void PCU_Comm_Order(bool on)
{
  if (global_state == uninit)
    pcu_fail("Comm_Order called before Comm_Init");
  pcu_msg* m = get_msg();
  if (on && (!m->order))
    m->order = pcu_order_new();
  if ((!on) && m->order) {
    pcu_order_free(m->order);
    m->order = NULL;
  }
}

/** \brief Blocking barrier over all threads. */
void PCU_Barrier(void)
{
  if (global_state == uninit)
    pcu_fail("Barrier called before Comm_Init");
  pcu_barrier(&(get_msg()->coll));
}

/** \brief Performs an Allreduce sum of double arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n doubles.
  After this call, p[i] will contain the sum of all p[i]'s
  given by each rank.
  */
void PCU_Add_Doubles(double* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Add_Doubles called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_doubles,p,n*sizeof(double));
}

/** \brief Performs an Allreduce minimum of double arrays.
  */
void PCU_Min_Doubles(double* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Min_Doubles called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_min_doubles,p,n*sizeof(double));
}

/** \brief Performs an Allreduce maximum of double arrays.
  */
void PCU_Max_Doubles(double* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Max_Doubles called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_max_doubles,p,n*sizeof(double));
}

/** \brief Performs an Allreduce sum of integers
  */
void PCU_Add_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Add_Ints called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_ints,p,n*sizeof(int));
}

/** \brief Performs an Allreduce sum of long integers
  */
void PCU_Add_Longs(long* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Add_Longs called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_longs,p,n*sizeof(long));
}

/** \brief Performs an exclusive prefix sum of integer arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n integers.
  After this call, p[i] will contain the sum of all p[i]'s
  given by ranks lower than the calling rank.
  */
void PCU_Exscan_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Exscan_Ints called before Comm_Init");
  int* originals;
  PCU_MALLOC(originals,n);
  for (size_t i=0; i < n; ++i)
    originals[i] = p[i];
  pcu_scan(&(get_msg()->coll),pcu_add_ints,p,n*sizeof(int));
  //convert inclusive scan to exclusive
  for (size_t i=0; i < n; ++i)
    p[i] -= originals[i];
  pcu_free(originals);
}

/** \brief See PCU_Exscan_Ints */
void PCU_Exscan_Longs(long* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Exscan_Longs called before Comm_Init");
  long* originals;
  PCU_MALLOC(originals,n);
  for (size_t i=0; i < n; ++i)
    originals[i] = p[i];
  pcu_scan(&(get_msg()->coll),pcu_add_longs,p,n*sizeof(long));
  //convert inclusive scan to exclusive
  for (size_t i=0; i < n; ++i)
    p[i] -= originals[i];
  pcu_free(originals);
}

/** \brief Performs an Allreduce minimum of int arrays.
  */
void PCU_Min_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Min_Ints called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_min_ints,p,n*sizeof(int));
}

/** \brief Performs an Allreduce maximum of int arrays.
  */
void PCU_Max_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    pcu_fail("Max_Ints called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_max_ints,p,n*sizeof(int));
}

/** \brief Performs a parallel logical OR reduction
  */
int PCU_Or(int c)
{
  PCU_Max_Ints(&c, 1);
  return c;
}

#if ENABLE_THREADS
static void* run(void* in)
{
  /* this wrapper around the user thread function
     sets up the PCU thread environment, including
     thread rank and thread-local messenger */
  pcu_thread_init(in);
  int rank = pcu_thread_rank();
  pcu_make_msg(global_tmsg + rank);
  if (global_args)
    global_args[rank] = global_function(global_args[rank]);
  else
    global_function(NULL);
  if (global_tmsg[rank].order)
    pcu_order_free(global_tmsg[rank].order);
  pcu_free_msg(global_tmsg + rank);
  return NULL;
}
#endif

/** \brief Runs \a nthreads instances of \a function, each in a thread.
  \details This function will create (\a nthreads - 1) new pthreads and use
  these as well as the caller thread to run \a function.
  The argument passed to thread i is \a in_out [i], and the return value
  of thread i is then stored in \a in_out [i].
  If in_out is NULL, all threads will receive NULL as their argument.

  Currently, PCU requires that this call is collective and homogeneous.
  This means that all processes in an MPI job should call PCU_Thrd_Run
  at the same time, and they should all pass the same number for \a nthreads.
  MPI_Init_thread should have been called before this function.

  Any calls to PCU_Comm functions from within one of these threads
  will have access to the hybrid communication interface.
  This means that ranks will be unique to a thread in the whole MPI job,
  and messages are sent and received between threads.
  Phases will be synchronized across all threads in the MPI job.
 */

int PCU_Thrd_Run(int nthreads, PCU_Thrd_Func function, void** in_out)
{
#if ENABLE_THREADS
  if (global_state == uninit)
    pcu_fail("Thrd_Run called before Comm_Init");
  if (pcu_get_mpi() == &pcu_tmpi)
    pcu_fail("nested calls to Thrd_Run");
  pcu_tmpi_check_support();
  if (global_pmsg.order)
    pcu_order_free(global_pmsg.order);
  pcu_free_msg(&global_pmsg);
  pcu_set_mpi(&pcu_tmpi);
  PCU_MALLOC(global_tmsg,(size_t)nthreads);
  global_function = function;
  global_args = in_out;
  pcu_run_threads(nthreads,run);
  pcu_free(global_tmsg);
  pcu_set_mpi(&pcu_pmpi);
  pcu_make_msg(&global_pmsg);
#else
  (void)nthreads;//unused parameter warning silencer
  (void)function;
  (void)in_out;
  fail_no_threads();
#endif
  return PCU_SUCCESS;
}

/** \brief Returns the process-unique rank of the calling thread.
  \details When called from a thread inside PCU_Thrd_Run, the resulting rank
  will be unique only within the same process.
  Ranks are contiguous integers from 0 to nthreads-1, with the thread that
  called PCU_Thrd_Run being assigned rank 0.
 */
int PCU_Thrd_Self(void)
{
#if ENABLE_THREADS
  if (pcu_get_mpi()==&pcu_tmpi)
    return pcu_thread_rank();
#endif
  return 0;
}

/** \brief Returns the number of threads running in the current process.
  \details When called from a thread inside PCU_Thrd_Run, returns the number
  of threads running in this process,
  which is equivalent to the nthreads argument to PCU_Thrd_Run.
 */
int PCU_Thrd_Peers(void)
{
#if ENABLE_THREADS
  if (pcu_get_mpi()==&pcu_tmpi)
    return pcu_thread_size();
#endif
  return 1;
}

/** \brief Blocks all threads of a process until all have hit the barrier.
 */
void PCU_Thrd_Barrier(void)
{
  if (global_state == uninit)
    pcu_fail("Thrd_Barrier called before Comm_Init");
#if ENABLE_THREADS  
  if (pcu_get_mpi()==&pcu_tmpi)
    pcu_thread_barrier();
#endif
}

/** \brief Acquire the PCU master spinlock.
  \details usage is discouraged, this function
  and PCU_Thrd_Unlock exist as a patch to allow
  people to use thread-unsafe third-party code
  at the expense of performance by serializing
  calls to said code.
 */
void PCU_Thrd_Lock(void)
{
#if ENABLE_THREADS  
  if (pcu_get_mpi()==&pcu_tmpi)
    pcu_thread_lock();
#endif
}

/** \brief Release the PCU master spinlock. */
void PCU_Thrd_Unlock(void)
{
#if ENABLE_THREADS  
  if (pcu_get_mpi()==&pcu_tmpi)
    pcu_thread_unlock();
#endif
}

/** \brief Returns the unique rank of the calling process.
 */
int PCU_Proc_Self(void)
{
  if (global_state == uninit)
    pcu_fail("Proc_Self called before Comm_Init");
  return pcu_pmpi_rank();
}

/** \brief Returns the number of processes.
 */
int PCU_Proc_Peers(void)
{
  if (global_state == uninit)
    pcu_fail("Proc_Peers called before Comm_Init");
  return pcu_pmpi_size();
}

/** \brief Similar to PCU_Comm_Self, returns the rank as an argument.
 */
int PCU_Comm_Rank(int* rank)
{
  if (global_state == uninit)
    pcu_fail("Comm_Rank called before Comm_Init");
  *rank = pcu_mpi_rank();
  return PCU_SUCCESS;
}

/** \brief Similar to PCU_Comm_Peers, returns the size as an argument. */
int PCU_Comm_Size(int* size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Size called before Comm_Init");
  *size = pcu_mpi_size();
  return PCU_SUCCESS;
}

/** \brief Returns true iff PCU has been initialized */
bool PCU_Comm_Initialized(void)
{
  return global_state == init;
}

/** \brief Deprecated, see PCU_Comm_Begin.
 */
int PCU_Comm_Start(PCU_Method method)
{
  (void)method; //warning silencer
  if (global_state == uninit)
    pcu_fail("Comm_Start called before Comm_Init");
  pcu_msg_start(get_msg());
  return PCU_SUCCESS;
}

/** \brief Returns in * \a size the number of bytes being sent to \a to_rank.
  \details Returns the size of the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Packed(int to_rank, size_t* size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Packed called before Comm_Init");
  if ((to_rank < 0)||(to_rank >= pcu_mpi_size()))
    pcu_fail("Invalid rank in Comm_Packed");
  *size = pcu_msg_packed(get_msg(),to_rank);
  return PCU_SUCCESS;
}

/** \brief Packs a message to be sent to \a to_rank.
  \details This function packs a message into the buffer being sent
  to \a to_rank.
  Messages packed by this function can be received using the function
  PCU_Comm_Read.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
  If this function is used, PCU_Comm_Pack should not be used.
 */
int PCU_Comm_Write(int to_rank, const void* data, size_t size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Write called before Comm_Init");
  if ((to_rank < 0)||(to_rank >= pcu_mpi_size()))
    pcu_fail("Invalid rank in Comm_Write");
  pcu_msg* msg = get_msg();
  PCU_MSG_PACK(msg,to_rank,size);
  memcpy(pcu_msg_pack(msg,to_rank,size),data,size);
  return PCU_SUCCESS;
}

/** \brief Convenience wrapper over Listen and Unpacked */
bool PCU_Comm_Receive(void)
{
  while (PCU_Comm_Unpacked())
    if (!PCU_Comm_Listen())
      return false;
  return true;
}

/** \brief Receives a message for this communication phase.
  \details This function tries to receive a message packed by
  PCU_Comm_Write.
  If a the communication phase is over and there are no more
  messages to receive, this function returns false.
  Otherwise, * \a from_rank will be the rank which sent the message,
 * \a data will point to the start of the message data, and
 * \a size will be the number of bytes of message data.
 If this function is used, PCU_Comm_Receive should not be used.
 Note that the address * \a data points into a PCU buffer, so
 it is strongly recommended that this data be read and not modified.
 */
bool PCU_Comm_Read(int* from_rank, void** data, size_t* size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Read called before Comm_Init");
  if (!PCU_Comm_Receive())
    return false;
  *from_rank = PCU_Comm_Sender();
  PCU_COMM_UNPACK(*size);
  *data = PCU_Comm_Extract(*size);
  return true;
}

/** \brief Open file debugN.txt, where N = PCU_Comm_Self(). */
void PCU_Debug_Open(void)
{
  if (global_state == uninit)
    pcu_fail("Debug_Open called before Comm_Init");
  pcu_msg* msg = get_msg();
  if ( ! msg->file)
    msg->file = pcu_open_parallel("debug","txt");
}

/** \brief like fprintf, contents go to debugN.txt */
void PCU_Debug_Print(const char* format, ...)
{
  if (global_state == uninit)
    pcu_fail("Debug_Print called before Comm_Init");
  pcu_msg* msg = get_msg();
  if ( ! msg->file)
    return; //Print is a no-op if no file is open
  va_list ap;
  va_start(ap,format);
  vfprintf(msg->file,format,ap);
  va_end(ap);
  fflush(msg->file);
}

/** \brief Similar to PCU_Comm_Sender, returns the rank as an argument. */
int PCU_Comm_From(int* from_rank)
{
  if (global_state == uninit)
    pcu_fail("Comm_From called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    *from_rank = pcu_order_received_from(m->order);
  else
    *from_rank = pcu_msg_received_from(m);
  return PCU_SUCCESS;
}

/** \brief Returns in * \a size the bytes in the current received buffer
  \details This function should be called after a successful PCU_Comm_Receive.
  The size returned will be the total received size regardless of how
  much unpacking has been done.
 */
int PCU_Comm_Received(size_t* size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Received called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    *size = pcu_order_received_size(m->order);
  else
    *size = pcu_msg_received_size(m);
  return PCU_SUCCESS;
}

/** \brief Extracts a block of data from the current received buffer.
  \details This function should be called after a successful PCU_Comm_Receive.
  The next \a size bytes of the current received buffer are unpacked,
  and an internal pointer to that data is returned.
  The returned pointer must not be freed by the user.
 */
void* PCU_Comm_Extract(size_t size)
{
  if (global_state == uninit)
    pcu_fail("Comm_Extract called before Comm_Init");
  pcu_msg* m = get_msg();
  if (m->order)
    return pcu_order_unpack(m->order,size);
  return pcu_msg_unpack(m,size);
}

/** \brief Reinitializes PCU with a new MPI communicator.
 \details All of PCU's logic is based off two duplicates
 of this communicator, so you can safely get PCU to act
 on sub-groups of processes using this function.
 This call should be collective over all processes
 in the previous communicator.
 Please do not mix this feature with PCU threading. 
 */
void PCU_Switch_Comm(MPI_Comm new_comm)
{
  if (global_state == uninit)
    pcu_fail("Switch_Comm called before Comm_Init");
  pcu_pmpi_switch(new_comm);
}

/** \brief Return the current MPI communicator
  \details Returns the communicator given to the 
  most recent PCU_Switch_Comm call, or MPI_COMM_WORLD
  otherwise.
 */
MPI_Comm PCU_Get_Comm(void)
{
  if (global_state == uninit)
    pcu_fail("Get_Comm called before Comm_Init");
  return pcu_pmpi_comm();
}

double PCU_Time(void)
{
  return MPI_Wtime();
}

