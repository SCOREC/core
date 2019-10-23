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
  based on MPI.
  PCU provides three things to users:
    1. A hybrid phased message passing system
    2. Hybrid collective operations

  Phased message passing is similar to Bulk Synchronous Parallel.
  All messages are exchanged in a phase, which is a collective operation
  involving all threads in the parallel program.
  During a phase, the following events happen in sequence:
    1. All threads send non-blocking messages to other threads
    2. All threads receive all messages sent to them during this phase
  PCU provides termination detection, which is the ability to detect when all
  messages have been received without prior knowledge of which threads
  are sending to which.

  The API documentation is here: pcu.c
*/

#include <string.h>
#include <stdarg.h>
#include "PCU.h"
#include "pcu_msg.h"
#include "pcu_pmpi.h"
#include "pcu_order.h"
#include "noto_malloc.h"
#include "reel.h"
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */

enum state { uninit, init };
static enum state global_state = uninit;
static pcu_msg global_pmsg;

static pcu_msg* get_msg()
{
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
    reel_fail("nested calls to Comm_Init");
  pcu_pmpi_init(MPI_COMM_WORLD);
  pcu_set_mpi(&pcu_pmpi);
  pcu_make_msg(&global_pmsg);
  global_state = init;
  /* turn ordering on by default, call
     PCU_Comm_Order(false) after PCU_Comm_Init
     to disable this */
  PCU_Comm_Order(true);
  return PCU_SUCCESS;
}

/** \brief Frees all PCU library structures.
  \details This function must be called by all MPI processes after all
  other calls to PCU, and before calling MPI_Finalize.
 */
int PCU_Comm_Free(void)
{
  if (global_state == uninit)
    reel_fail("Comm_Free called before Comm_Init");
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

  Ranks are consecutive from 0 to \f$pt-1\f$ for a program with
  \f$p\f$ processes and \f$t\f$ threads per process.
  Ranks are contiguous within a process, so that the \f$t\f$ threads in process
  \f$i\f$ are numbered from \f$ti\f$ to \f$ti+t-1\f$.
 */
int PCU_Comm_Self(void)
{
  if (global_state == uninit)
    reel_fail("Comm_Self called before Comm_Init");
  return pcu_mpi_rank();
}

/** \brief Returns the number of threads in the program.
  \details when called from a non-threaded MPI process, this function is
  equivalent to MPI_Comm_size(MPI_COMM_WORLD,size).
 */
int PCU_Comm_Peers(void)
{
  if (global_state == uninit)
    reel_fail("Comm_Peers called before Comm_Init");
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
    reel_fail("Comm_Begin called before Comm_Init");
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
    reel_fail("Comm_Pack called before Comm_Init");
  if ((to_rank < 0)||(to_rank >= pcu_mpi_size()))
    reel_fail("Invalid rank in Comm_Pack");
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
    reel_fail("Comm_Send called before Comm_Init");
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
    reel_fail("Comm_Listen called before Comm_Init");
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
    reel_fail("Comm_Sender called before Comm_Init");
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
    reel_fail("Comm_Unpacked called before Comm_Init");
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
    reel_fail("Comm_Unpack called before Comm_Init");
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
    reel_fail("Comm_Order called before Comm_Init");
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
    reel_fail("Barrier called before Comm_Init");
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
    reel_fail("Add_Doubles called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_doubles,p,n*sizeof(double));
}

double PCU_Add_Double(double x)
{
  double a[1];
  a[0] = x;
  PCU_Add_Doubles(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce minimum of double arrays.
  */
void PCU_Min_Doubles(double* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Min_Doubles called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_min_doubles,p,n*sizeof(double));
}

double PCU_Min_Double(double x)
{
  double a[1];
  a[0] = x;
  PCU_Min_Doubles(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce maximum of double arrays.
  */
void PCU_Max_Doubles(double* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Max_Doubles called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_max_doubles,p,n*sizeof(double));
}

double PCU_Max_Double(double x)
{
  double a[1];
  a[0] = x;
  PCU_Max_Doubles(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce sum of integers
  */
void PCU_Add_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Add_Ints called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_ints,p,n*sizeof(int));
}

int PCU_Add_Int(int x)
{
  int a[1];
  a[0] = x;
  PCU_Add_Ints(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce sum of long integers
  */
void PCU_Add_Longs(long* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Add_Longs called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_longs,p,n*sizeof(long));
}

long PCU_Add_Long(long x)
{
  long a[1];
  a[0] = x;
  PCU_Add_Longs(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce sum of size_t unsigned integers
  */
void PCU_Add_SizeTs(size_t* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Add_SizeTs called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_add_sizets,p,n*sizeof(size_t));
}

size_t PCU_Add_SizeT(size_t x)
{
  size_t a[1];
  a[0] = x;
  PCU_Add_SizeTs(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce minimum of size_t unsigned integers
  */
void PCU_Min_SizeTs(size_t* p, size_t n) {
  if (global_state == uninit)
    reel_fail("Min_SizeTs called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_min_sizets,p,n*sizeof(size_t));
}

size_t PCU_Min_SizeT(size_t x) {
  size_t a[1];
  a[0] = x;
  PCU_Min_SizeTs(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce maximum of size_t unsigned integers
  */
void PCU_Max_SizeTs(size_t* p, size_t n) {
  if (global_state == uninit)
    reel_fail("Max_SizeTs called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_max_sizets,p,n*sizeof(size_t));
}

size_t PCU_Max_SizeT(size_t x) {
  size_t a[1];
  a[0] = x;
  PCU_Max_SizeTs(a, 1);
  return a[0];
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
    reel_fail("Exscan_Ints called before Comm_Init");
  int* originals;
  NOTO_MALLOC(originals,n);
  for (size_t i=0; i < n; ++i)
    originals[i] = p[i];
  pcu_scan(&(get_msg()->coll),pcu_add_ints,p,n*sizeof(int));
  //convert inclusive scan to exclusive
  for (size_t i=0; i < n; ++i)
    p[i] -= originals[i];
  noto_free(originals);
}

int PCU_Exscan_Int(int x)
{
  int a[1];
  a[0] = x;
  PCU_Exscan_Ints(a, 1);
  return a[0];
}

/** \brief See PCU_Exscan_Ints */
void PCU_Exscan_Longs(long* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Exscan_Longs called before Comm_Init");
  long* originals;
  NOTO_MALLOC(originals,n);
  for (size_t i=0; i < n; ++i)
    originals[i] = p[i];
  pcu_scan(&(get_msg()->coll),pcu_add_longs,p,n*sizeof(long));
  //convert inclusive scan to exclusive
  for (size_t i=0; i < n; ++i)
    p[i] -= originals[i];
  noto_free(originals);
}

long PCU_Exscan_Long(long x)
{
  long a[1];
  a[0] = x;
  PCU_Exscan_Longs(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce minimum of int arrays.
  */
void PCU_Min_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Min_Ints called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_min_ints,p,n*sizeof(int));
}

int PCU_Min_Int(int x)
{
  int a[1];
  a[0] = x;
  PCU_Min_Ints(a, 1);
  return a[0];
}

/** \brief Performs an Allreduce maximum of int arrays.
  */
void PCU_Max_Ints(int* p, size_t n)
{
  if (global_state == uninit)
    reel_fail("Max_Ints called before Comm_Init");
  pcu_allreduce(&(get_msg()->coll),pcu_max_ints,p,n*sizeof(int));
}

int PCU_Max_Int(int x)
{
  int a[1];
  a[0] = x;
  PCU_Max_Ints(a, 1);
  return a[0];
}

/** \brief Performs a parallel logical OR reduction
  */
int PCU_Or(int c)
{
  return PCU_Max_Int(c);
}

/** \brief Performs a parallel logical AND reduction
  */
int PCU_And(int c)
{
  return PCU_Min_Int(c);
}

/** \brief Performs a parallel gather followed by a broadcast
  */
void PCU_Allgather_Doubles(double* i, size_t ni, double* o)
{
  if(global_state == uninit)
    reel_fail("Allgather_Doubles called before Comm_Init");
  MPI_Allgather(i,ni,MPI_DOUBLE,o,ni,MPI_DOUBLE,PCU_Get_Comm());
}

void PCU_Allgather_Double(double i, double* o)
{
  PCU_Allgather_Doubles(&i,1,o);
}

/** \brief Performs a parallel gather followed by a broadcast
  */
void PCU_Allgather_Ints(int* i, size_t ni, int* o)
{
  if(global_state == uninit)
    reel_fail("Allgather_Ints called before Comm_Init");
  MPI_Allgather(i,ni,MPI_INTEGER,o,ni,MPI_INTEGER,PCU_Get_Comm());
}

void PCU_Allgather_Int(int i, int* o)
{
  PCU_Allgather_Ints(&i,1,o);
}

/** \brief Performs a parallel gather followed by a broadcast
  */
void PCU_Allgather_Longs(long* i, size_t ni, long* o)
{
  if(global_state == uninit)
    reel_fail("Allgather_Longs called before Comm_Init");
  MPI_Allgather(i,ni,MPI_LONG,o,ni,MPI_LONG,PCU_Get_Comm());
}

void PCU_Allgather_Long(long i, long* o)
{
  PCU_Allgather_Longs(&i,1,o);
}

/** \brief Performs a parallel gather followed by a broadcast
  */
void PCU_Allgather_SizeTs(size_t* i, size_t ni, size_t *o)
{
  if(global_state == uninit)
    reel_fail("Allgather_SizeTs called before Comm_Init");
  MPI_Allgather(i,ni*sizeof(size_t),MPI_BYTE,o,ni*sizeof(size_t),MPI_BYTE,PCU_Get_Comm());
}

void PCU_Allgather_SizeT(size_t i, size_t* o)
{
  PCU_Allgather_SizeTs(&i,1,o);
}

/** \brief Returns the unique rank of the calling process.
 */
int PCU_Proc_Self(void)
{
  if (global_state == uninit)
    reel_fail("Proc_Self called before Comm_Init");
  return pcu_pmpi_rank();
}

/** \brief Returns the number of processes.
 */
int PCU_Proc_Peers(void)
{
  if (global_state == uninit)
    reel_fail("Proc_Peers called before Comm_Init");
  return pcu_pmpi_size();
}

/** \brief Similar to PCU_Comm_Self, returns the rank as an argument.
 */
int PCU_Comm_Rank(int* rank)
{
  if (global_state == uninit)
    reel_fail("Comm_Rank called before Comm_Init");
  *rank = pcu_mpi_rank();
  return PCU_SUCCESS;
}

/** \brief Similar to PCU_Comm_Peers, returns the size as an argument. */
int PCU_Comm_Size(int* size)
{
  if (global_state == uninit)
    reel_fail("Comm_Size called before Comm_Init");
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
    reel_fail("Comm_Start called before Comm_Init");
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
    reel_fail("Comm_Packed called before Comm_Init");
  if ((to_rank < 0)||(to_rank >= pcu_mpi_size()))
    reel_fail("Invalid rank in Comm_Packed");
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
    reel_fail("Comm_Write called before Comm_Init");
  if ((to_rank < 0)||(to_rank >= pcu_mpi_size()))
    reel_fail("Invalid rank in Comm_Write");
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
    reel_fail("Comm_Read called before Comm_Init");
  if (!PCU_Comm_Receive())
    return false;
  *from_rank = PCU_Comm_Sender();
  PCU_COMM_UNPACK(*size);
  *data = PCU_Comm_Extract(*size);
  return true;
}

static void safe_mkdir(const char* path, mode_t mode)
{
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
    reel_fail("PCU: could not create directory \"%s\"\n", path);
}

static void append(char* s, size_t size, const char* format, ...)
{
  int len = strlen(s);
  va_list ap;
  va_start(ap, format);
  vsnprintf(s + len, size - len, format, ap);
  va_end(ap);
}


void PCU_Debug_Open(void)
{
  if (global_state == uninit)
    reel_fail("Debug_Open called before Comm_Init");

  const int fanout = 2048;
  const int bufsize = 1024;
  char* path = noto_malloc(bufsize);
  path[0] = '\0';
  if (PCU_Comm_Peers() > fanout) {
    mode_t const dir_perm = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
    strcpy(path, "debug/");
    safe_mkdir(path, dir_perm);
    int self = PCU_Comm_Self();
    append(path, bufsize, "%d/", self / fanout);
    if (self % fanout == 0)
      safe_mkdir(path, dir_perm);
    PCU_Barrier();
  }

  append(path,bufsize, "%s", "debug");
  pcu_msg* msg = get_msg();
  if ( ! msg->file)
    msg->file = pcu_open_parallel(path,"txt");
  noto_free(path);
}

/** \brief like fprintf, contents go to debugN.txt */
void PCU_Debug_Print(const char* format, ...)
{
  if (global_state == uninit)
    reel_fail("Debug_Print called before Comm_Init");
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
    reel_fail("Comm_From called before Comm_Init");
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
    reel_fail("Comm_Received called before Comm_Init");
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
    reel_fail("Comm_Extract called before Comm_Init");
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
 */
void PCU_Switch_Comm(MPI_Comm new_comm)
{
  if (global_state == uninit)
    reel_fail("Switch_Comm called before Comm_Init");
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
    reel_fail("Get_Comm called before Comm_Init");
  return pcu_pmpi_comm();
}

int PCU_Local_To_Foreign(int local_rank, MPI_Comm foreign_comm)
{
  if (global_state == uninit)
    reel_fail("Local_To_Foreign called before Comm_Init");
  return pcu_pmpi_lcl_to_frn(local_rank,foreign_comm);
}

int PCU_Foreign_To_Local(int foreign_rank, MPI_Comm foreign_comm)
{
  if (global_state == uninit)
    reel_fail("Foreign_To_Local called before Comm_Init");
  return pcu_pmpi_frn_to_lcl(foreign_rank,foreign_comm);
}

/** \brief Return the time in seconds since some time in the past
 */
double PCU_Time(void)
{
  return MPI_Wtime();
}

void PCU_Protect(void)
{
  reel_protect();
}
