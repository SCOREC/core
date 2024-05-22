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
#include "PCUObj.h"
#include "pcu_msg.h"
#if defined(SCOREC_NO_MPI)
  #include "pcu_pnompi.h"
#else
  #include "pcu_pmpi.h"
#endif
#include "pcu_order.h"
#include "reel.h"
#include <cstdarg>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
#include <limits.h> /*INT_MAX*/
#include <stdlib.h> /*abort*/

static pcu::PCU *global_pcu = nullptr;
namespace pcu {
  pcu::PCU* PCU_GetGlobal() { return global_pcu; }
}

extern "C" {

/** \brief Initializes the PCU library.
  \details This function must be called by all MPI processes before
  calling any other PCU functions.
  MPI_Init or MPI_Init_thread should be called before this function.
 */
int PCU_Comm_Init(void) {
  if (global_pcu != nullptr)
    reel_fail("nested calls to Comm_Init");
  global_pcu = new pcu::PCU(MPI_COMM_WORLD);
  return PCU_SUCCESS;
}

/** \brief Frees all PCU library structures.
  \details This function must be called by all MPI processes after all
  other calls to PCU, and before calling MPI_Finalize.
 */
int PCU_Comm_Free(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Free called before Comm_Init");
  delete global_pcu;
  global_pcu = nullptr;
  return PCU_SUCCESS;
}

int PCU_Comm_Free_One(MPI_Comm* com)
{
  if (global_pcu == nullptr)
    reel_fail("Comm_Free_One called before Comm_Init");
  return global_pcu->Free_One(com);
}

int PCU_Comm_Split(MPI_Comm oldCom, int color, int key, MPI_Comm* newCom)
{
  if (global_pcu == nullptr)
    reel_fail("Comm_Split called before Comm_Init");
  return global_pcu->Split(oldcom, color, key, newCom);
}

int PCU_Comm_Allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (global_pcu == nullptr)
    reel_fail("Comm_Allreduce called before Comm_Init");
  return global_pcu->Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  pcu_pmpi_allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  return PCU_SUCCESS;
}

int PCU_Comm_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  pcu_pmpi_allgather(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm);
  return PCU_SUCCESS;
}

 int PCU_Comm_Barrier(MPI_Comm comm)
{
  pcu_pmpi_barrier(comm);
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
int PCU_Comm_Self(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Self called before Comm_Init");
  return global_pcu->Self();
}

/** \brief Returns the number of threads in the program.
  \details when called from a non-threaded MPI process, this function is
  equivalent to MPI_Comm_size(MPI_COMM_WORLD,size).
 */
int PCU_Comm_Peers(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Peers called before Comm_Init");
  return global_pcu->Peers();
}

/** \brief Begins a PCU communication phase.
  \details This function must be called by all threads in the MPI job
  at the beginning of each phase of communication.
  After calling this function, each thread may call functions like
  PCU_Comm_Pack or PCU_Comm_Write.
*/
void PCU_Comm_Begin(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Begin called before Comm_Init");
  global_pcu->Begin();
}

/** \brief Packs data to be sent to \a to_rank.
  \details This function appends the block of \a size bytes starting
  at \a data to the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Pack(int to_rank, const void *data, size_t size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Pack called before Comm_Init");
  return global_pcu->Pack(to_rank, data, size);
}

/** \brief Sends all buffers for this communication phase.
  \details This function should be called by all threads in the MPI job
  after calls to PCU_Comm_Pack or PCU_Comm_Write and before calls
  to PCU_Comm_Listen or PCU_Comm_Read.
  All buffers from this thread are sent out and receiving
  may begin after this call.
 */
int PCU_Comm_Send(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Send called before Comm_Init");
  return global_pcu->Send();
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
bool PCU_Comm_Listen(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Listen called before Comm_Init");
  return global_pcu->Listen();
}

/** \brief Returns in * \a from_rank the sender of the current received buffer.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
int PCU_Comm_Sender(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Sender called before Comm_Init");
  return global_pcu->Sender();
}

/** \brief Returns true if the current received buffer has been unpacked.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
bool PCU_Comm_Unpacked(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Unpacked called before Comm_Init");
  return global_pcu->Unpacked();
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
int PCU_Comm_Unpack(void *data, size_t size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Unpack called before Comm_Init");
  return global_pcu->Unpack(data, size);
}

void PCU_Comm_Order(bool on) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Order called before Comm_Init");
  global_pcu->Order(on);
}

/** \brief Blocking barrier over all threads. */
void PCU_Barrier(void) {
  if (global_pcu == nullptr)
    reel_fail("Barrier called before Comm_Init");
  global_pcu->Barrier();
}

/** \brief Performs an Allreduce sum of double arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n doubles.
  After this call, p[i] will contain the sum of all p[i]'s
  given by each rank.
  */
void PCU_Add_Doubles(double *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Add_Doubles called before Comm_Init");
  global_pcu->Add(p, n);
}

double PCU_Add_Double(double x) {
  if (global_pcu == nullptr)
    reel_fail("Add_Double called before Comm_Init");
  return global_pcu->Add(x);
}

/** \brief Performs an Allreduce minimum of double arrays.
 */
void PCU_Min_Doubles(double *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Min_Doubles called before Comm_Init");
  global_pcu->Min(p, n);
}

double PCU_Min_Double(double x) {
  if (global_pcu == nullptr)
    reel_fail("Min_Double called before Comm_Init");
  return global_pcu->Min(x);
}

/** \brief Performs an Allreduce maximum of double arrays.
 */
void PCU_Max_Doubles(double *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Max_Doubles called before Comm_Init");
  global_pcu->Max(p, n);
}

double PCU_Max_Double(double x) {
  if (global_pcu == nullptr)
    reel_fail("Max_Double called before Comm_Init");
  return global_pcu->Max(x);
}

/** \brief Performs an Allreduce sum of integers
 */
void PCU_Add_Ints(int *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Add_Ints called before Comm_Init");
  global_pcu->Add(p, n);
}

int PCU_Add_Int(int x) {
  if (global_pcu == nullptr)
    reel_fail("Add_Int called before Comm_Init");
  return global_pcu->Add(x);
}

/** \brief Performs an Allreduce sum of long integers
 */
void PCU_Add_Longs(long *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Add_Longs called before Comm_Init");
  global_pcu->Add(p, n);
}

long PCU_Add_Long(long x) {
  if (global_pcu == nullptr)
    reel_fail("Add_Long called before Comm_Init");
  return global_pcu->Add(x);
}

/** \brief Performs an Allreduce sum of size_t unsigned integers
 */
void PCU_Add_SizeTs(size_t *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Add_SizeTs called before Comm_Init");
  global_pcu->Add(p, n);
}

size_t PCU_Add_SizeT(size_t x) {
  if (global_pcu == nullptr)
    reel_fail("Add_SizeT called before Comm_Init");
  return global_pcu->Add(x);
}

/** \brief Performs an Allreduce minimum of size_t unsigned integers
 */
void PCU_Min_SizeTs(size_t *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Min_SizeTs called before Comm_Init");
  global_pcu->Min(p, n);
}

size_t PCU_Min_SizeT(size_t x) {
  if (global_pcu == nullptr)
    reel_fail("Min_SizeT called before Comm_Init");
  return global_pcu->Min(x);
}

/** \brief Performs an Allreduce maximum of size_t unsigned integers
 */
void PCU_Max_SizeTs(size_t *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Max_SizeTs called before Comm_Init");
  global_pcu->Max(p, n);
}

size_t PCU_Max_SizeT(size_t x) {
  if (global_pcu == nullptr)
    reel_fail("Max_SizeT called before Comm_Init");
  return global_pcu->Max(x);
}

/** \brief Performs an exclusive prefix sum of integer arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n integers.
  After this call, p[i] will contain the sum of all p[i]'s
  given by ranks lower than the calling rank.
  */
void PCU_Exscan_Ints(int *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Exscan_Ints called before Comm_Init");
  global_pcu->Exscan(p, n);
}

int PCU_Exscan_Int(int x) {
  if (global_pcu == nullptr)
    reel_fail("Exscan_Int called before Comm_Init");
  return global_pcu->Exscan(x);
}

/** \brief See PCU_Exscan_Ints */
void PCU_Exscan_Longs(long *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Exscan_Longs called before Comm_Init");
  global_pcu->Exscan(p, n);
}

long PCU_Exscan_Long(long x) {
  if (global_pcu == nullptr)
    reel_fail("Exscan_Long called before Comm_Init");
  return global_pcu->Exscan(x);
}

/** \brief Performs an Allreduce minimum of int arrays.
 */
void PCU_Min_Ints(int *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Min_Ints called before Comm_Init");
  global_pcu->Min(p, n);
}

int PCU_Min_Int(int x) {
  if (global_pcu == nullptr)
    reel_fail("Min_Int called before Comm_Init");
  return global_pcu->Min(x);
}

/** \brief Performs an Allreduce maximum of int arrays.
 */
void PCU_Max_Ints(int *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Max_Ints called before Comm_Init");
  global_pcu->Max(p, n);
}

int PCU_Max_Int(int x) {
  if (global_pcu == nullptr)
    reel_fail("Max_Int called before Comm_Init");
  return global_pcu->Max(x);
}
/** \brief Performs an Allreduce maximum of long arrays.
 */
void PCU_Max_Longs(long *p, size_t n) {
  if (global_pcu == nullptr)
    reel_fail("Max_Longs called before Comm_Init");
  global_pcu->Max(p, n);
}

long PCU_Max_Long(long x) {
  if (global_pcu == nullptr)
    reel_fail("Max_Long called before Comm_Init");
  return global_pcu->Max(x);
}

/** \brief Performs a parallel logical OR reduction
 */
int PCU_Or(int c) {
  if (global_pcu == nullptr)
    reel_fail("Or called before Comm_Init");
  return global_pcu->Or(c);
}

/** \brief Performs a parallel logical AND reduction
 */
int PCU_And(int c) {
  if (global_pcu == nullptr)
    reel_fail("And called before Comm_Init");
  return global_pcu->And(c);
}

/** \brief Returns the unique rank of the calling process.
 */
int PCU_Proc_Self(void) {
  if (global_pcu == nullptr)
    reel_fail("Proc_Self called before Comm_Init");
  return global_pcu->Self();
}

/** \brief Returns the number of processes.
 */
int PCU_Proc_Peers(void) {
  if (global_pcu == nullptr)
    reel_fail("Proc_Peers called before Comm_Init");
  return global_pcu->Peers();
}

/** \brief Similar to PCU_Comm_Self, returns the rank as an argument.
 */
int PCU_Comm_Rank(int *rank) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Rank called before Comm_Init");
  *rank = global_pcu->Self();
  return PCU_SUCCESS;
}

/** \brief Similar to PCU_Comm_Peers, returns the size as an argument. */
int PCU_Comm_Size(int *size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Size called before Comm_Init");
  *size = global_pcu->Peers();
  return PCU_SUCCESS;
}

/** \brief Returns true iff PCU has been initialized */
bool PCU_Comm_Initialized(void) { return global_pcu != nullptr; }

/** \brief Deprecated, see PCU_Comm_Begin.
 */
int PCU_Comm_Start(PCU_Method method) {
  (void)method; // warning silencer
  global_pcu->Begin();
  return PCU_SUCCESS;
}

/** \brief Returns in * \a size the number of bytes being sent to \a to_rank.
  \details Returns the size of the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Packed(int to_rank, size_t *size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Packed called before Comm_Init");
  return global_pcu->Packed(to_rank, size);
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
int PCU_Comm_Write(int to_rank, const void *data, size_t size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Write called before Comm_Init");
  return global_pcu->Write(to_rank, data, size);
}

/** \brief Convenience wrapper over Listen and Unpacked */
bool PCU_Comm_Receive(void) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Receive called before Comm_Init");
  return global_pcu->Receive();
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
bool PCU_Comm_Read(int *from_rank, void **data, size_t *size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Read called before Comm_Init");
  return global_pcu->Read(from_rank, data, size);
}

void PCU_Debug_Open(void) {
  if (global_pcu == nullptr)
    reel_fail("Debug_Open called before Comm_Init");
  global_pcu->DebugOpen();
}

/** \brief like fprintf, contents go to debugN.txt */
void PCU_Debug_Print(const char *format, ...) {
  if (global_pcu == nullptr)
    reel_fail("Debug_Print called before Comm_Init");
  va_list arglist;
  va_start(arglist, format);
  global_pcu->DebugPrint(format, arglist);
  va_end(arglist);
}

/** \brief Similar to PCU_Comm_Sender, returns the rank as an argument. */
int PCU_Comm_From(int *from_rank) {
  if (global_pcu == nullptr)
    reel_fail("Comm_From called before Comm_Init");
  return global_pcu->From(from_rank);
}

/** \brief Returns in * \a size the bytes in the current received buffer
  \details This function should be called after a successful PCU_Comm_Receive.
  The size returned will be the total received size regardless of how
  much unpacking has been done.
 */
int PCU_Comm_Received(size_t *size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Received called before Comm_Init");
  return global_pcu->Received(size);
}

/** \brief Extracts a block of data from the current received buffer.
  \details This function should be called after a successful PCU_Comm_Receive.
  The next \a size bytes of the current received buffer are unpacked,
  and an internal pointer to that data is returned.
  The returned pointer must not be freed by the user.
 */
void *PCU_Comm_Extract(size_t size) {
  if (global_pcu == nullptr)
    reel_fail("Comm_Extract called before Comm_Init");
  return global_pcu->Extract(size);
}

/** \brief Reinitializes PCU with a new MPI communicator.
 \details All of PCU's logic is based off two duplicates
 of this communicator, so you can safely get PCU to act
 on sub-groups of processes using this function.
 This call should be collective over all processes
 in the previous communicator. This is a very heavy weight function
 and should be used sparingly.
 */
void PCU_Switch_Comm(MPI_Comm new_comm) {
  if (global_pcu == nullptr)
    reel_fail("Switch_Comm called before Comm_Init");
  global_pcu->SwitchMPIComm(new_comm);
}

/** \brief Return the current MPI communicator
  \details Returns the communicator given to the
  most recent PCU_Switch_Comm call, or MPI_COMM_WORLD
  otherwise.
 */
MPI_Comm PCU_Get_Comm(void) {
  if (global_pcu == nullptr)
    reel_fail("Get_Comm called before Comm_Init");
  return global_pcu->GetMPIComm();
}

/** \brief Return the time in seconds since some time in the past
 */
double PCU_Time(void) { return pcu::Time(); }

void PCU_Protect(void) { return pcu::Protect(); }

double PCU_GetMem(void) { return pcu::GetMem(); }


/** \brief Return the global PCU
 */

PCUHandle PCU_Get_Global_Handle(void) { 
  PCUHandle h;
  h.ptr = static_cast<void*>(pcu::PCU_GetGlobal());
  return h;
}

}
