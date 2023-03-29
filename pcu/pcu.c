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
#include "PCU2.h"
#include "pcu_msg.h"
#include "pcu_order.h"
#include "noto_malloc.h"
#include "reel.h"
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
#include <limits.h> /*INT_MAX*/
#include <stdlib.h> /*abort*/


static pcu_t global_pcu = {.state=pcu_state_uninit};

/** \brief Initializes the PCU library.
  \details This function must be called by all MPI processes before
  calling any other PCU functions.
  MPI_Init or MPI_Init_thread should be called before this function.
 */
int PCU_Comm_Init(void)
{
  PCU_Comm_Init_2(&global_pcu);
}

/** \brief Frees all PCU library structures.
  \details This function must be called by all MPI processes after all
  other calls to PCU, and before calling MPI_Finalize.
 */
int PCU_Comm_Free(void)
{
  return PCU_Comm_Free_2(&global_pcu);
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
  return PCU_Comm_Self_2(&global_pcu);
}

/** \brief Returns the number of threads in the program.
  \details when called from a non-threaded MPI process, this function is
  equivalent to MPI_Comm_size(MPI_COMM_WORLD,size).
 */
int PCU_Comm_Peers(void)
{
  return PCU_Comm_Peers_2(&global_pcu);
}

/** \brief Begins a PCU communication phase.
  \details This function must be called by all threads in the MPI job
  at the beginning of each phase of communication.
  After calling this function, each thread may call functions like
  PCU_Comm_Pack or PCU_Comm_Write.
*/
void PCU_Comm_Begin(void)
{
  PCU_Comm_Begin_2(&global_pcu);
}

/** \brief Packs data to be sent to \a to_rank.
  \details This function appends the block of \a size bytes starting
  at \a data to the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Pack(int to_rank, const void* data, size_t size)
{
  return PCU_Comm_Pack_2(&global_pcu, to_rank, data, size);
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
  return PCU_Comm_Send_2(&global_pcu);
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
  return PCU_Comm_Listen_2(&global_pcu);
}

/** \brief Returns in * \a from_rank the sender of the current received buffer.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
int PCU_Comm_Sender(void)
{
  return PCU_Comm_Sender_2(&global_pcu);
}

/** \brief Returns true if the current received buffer has been unpacked.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
bool PCU_Comm_Unpacked(void)
{
  return PCU_Comm_Unpacked_2(&global_pcu);
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
  return PCU_Comm_Unpack_2(&global_pcu, data, size);
}

void PCU_Comm_Order(bool on)
{
  PCU_Comm_Order_2(&global_pcu, on);
}

/** \brief Blocking barrier over all threads. */
void PCU_Barrier(void)
{
  PCU_Barrier_2(&global_pcu);
}

/** \brief Performs an Allreduce sum of double arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n doubles.
  After this call, p[i] will contain the sum of all p[i]'s
  given by each rank.
  */
void PCU_Add_Doubles(double* p, size_t n)
{
  PCU_Add_Doubles_2(&global_pcu, p, n);
}

double PCU_Add_Double(double x)
{
  return PCU_Add_Double_2(&global_pcu, x);
}

/** \brief Performs an Allreduce minimum of double arrays.
  */
void PCU_Min_Doubles(double* p, size_t n)
{
  PCU_Min_Doubles_2(&global_pcu, p, n);
}

double PCU_Min_Double(double x)
{
  return PCU_Min_Double_2(&global_pcu, x);
}

/** \brief Performs an Allreduce maximum of double arrays.
  */
void PCU_Max_Doubles(double* p, size_t n)
{
  PCU_Max_Doubles_2(&global_pcu, p, n);
}

double PCU_Max_Double(double x)
{
  return PCU_Max_Double_2(&global_pcu, x);
}

/** \brief Performs an Allreduce sum of integers
  */
void PCU_Add_Ints(int* p, size_t n)
{
  PCU_Add_Ints_2(&global_pcu, p, n);
}

int PCU_Add_Int(int x)
{
  return PCU_Add_Int_2(&global_pcu, x);
}

/** \brief Performs an Allreduce sum of long integers
  */
void PCU_Add_Longs(long* p, size_t n)
{
  PCU_Add_Longs_2(&global_pcu, p, n);
}

long PCU_Add_Long(long x)
{
  return PCU_Add_Long_2(&global_pcu, x);
}

/** \brief Performs an Allreduce sum of size_t unsigned integers
  */
void PCU_Add_SizeTs(size_t* p, size_t n)
{
  PCU_Add_SizeTs_2(&global_pcu, p, n);
}

size_t PCU_Add_SizeT(size_t x)
{
  return PCU_Add_SizeT_2(&global_pcu, x);
}

/** \brief Performs an Allreduce minimum of size_t unsigned integers
  */
void PCU_Min_SizeTs(size_t* p, size_t n) {
  PCU_Min_SizeTs_2(&global_pcu, p, n);
}

size_t PCU_Min_SizeT(size_t x) {
  return PCU_Min_SizeT_2(&global_pcu, x);
}

/** \brief Performs an Allreduce maximum of size_t unsigned integers
  */
void PCU_Max_SizeTs(size_t* p, size_t n) {
  PCU_Max_SizeTs_2(&global_pcu, p, n);
}

size_t PCU_Max_SizeT(size_t x) {
  return PCU_Max_SizeT_2(&global_pcu, x);
}

/** \brief Performs an exclusive prefix sum of integer arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n integers.
  After this call, p[i] will contain the sum of all p[i]'s
  given by ranks lower than the calling rank.
  */
void PCU_Exscan_Ints(int* p, size_t n)
{
  PCU_Exscan_Ints_2(&global_pcu, p, n);
}

int PCU_Exscan_Int(int x)
{
  return PCU_Exscan_Int_2(&global_pcu, x);
}

/** \brief See PCU_Exscan_Ints */
void PCU_Exscan_Longs(long* p, size_t n)
{
  PCU_Exscan_Longs_2(&global_pcu, p, n);
}

long PCU_Exscan_Long(long x)
{
  return PCU_Exscan_Long_2(&global_pcu, x);
}

/** \brief Performs an Allreduce minimum of int arrays.
  */
void PCU_Min_Ints(int* p, size_t n)
{
  PCU_Min_Ints_2(&global_pcu, p, n);
}

int PCU_Min_Int(int x)
{
  return PCU_Min_Int_2(&global_pcu, x);
}

/** \brief Performs an Allreduce maximum of int arrays.
  */
void PCU_Max_Ints(int* p, size_t n)
{
  PCU_Max_Ints_2(&global_pcu, p, n);
}

int PCU_Max_Int(int x)
{
  return PCU_Max_Int_2(&global_pcu, x);
}
/** \brief Performs an Allreduce maximum of long arrays.
  */
void PCU_Max_Longs(long* p, size_t n)
{
  PCU_Max_Longs_2(&global_pcu, p, n);
}

long PCU_Max_Long(long x)
{
  return PCU_Max_Long_2(&global_pcu, x);
}

/** \brief Performs a parallel logical OR reduction
  */
int PCU_Or(int c)
{
  return PCU_Or_2(&global_pcu, c);
}

/** \brief Performs a parallel logical AND reduction
  */
int PCU_And(int c)
{
  return PCU_And_2(&global_pcu, c);
}

/** \brief Returns the unique rank of the calling process.
 */
int PCU_Proc_Self(void)
{
  return PCU_Proc_Self_2(&global_pcu);
}

/** \brief Returns the number of processes.
 */
int PCU_Proc_Peers(void)
{
  return PCU_Proc_Peers_2(&global_pcu);
}

/** \brief Similar to PCU_Comm_Self, returns the rank as an argument.
 */
int PCU_Comm_Rank(int* rank)
{
  return PCU_Comm_Rank_2(&global_pcu, rank);
}

/** \brief Similar to PCU_Comm_Peers, returns the size as an argument. */
int PCU_Comm_Size(int* size)
{
  return PCU_Comm_Size_2(&global_pcu, size);
}

/** \brief Returns true iff PCU has been initialized */
bool PCU_Comm_Initialized(void)
{
  return PCU_Comm_Initialized_2(&global_pcu);
}

/** \brief Deprecated, see PCU_Comm_Begin.
 */
int PCU_Comm_Start(PCU_Method method)
{
  (void)method; //warning silencer
  return PCU_Comm_Start_2(&global_pcu);
}

/** \brief Returns in * \a size the number of bytes being sent to \a to_rank.
  \details Returns the size of the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Packed(int to_rank, size_t* size)
{
  return PCU_Comm_Packed_2(&global_pcu, to_rank, size);
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
  return PCU_Comm_Write_2(&global_pcu, to_rank, data, size);
}

/** \brief Convenience wrapper over Listen and Unpacked */
bool PCU_Comm_Receive(void)
{
  return PCU_Comm_Receive_2(&global_pcu);
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
  return PCU_Comm_Read_2(&global_pcu, from_rank, data, size);
}




void PCU_Debug_Open(void)
{
  PCU_Debug_Open_2(&global_pcu);
}

/** \brief like fprintf, contents go to debugN.txt */
void PCU_Debug_Print(const char* format, ...)
{
  va_list arglist;
  va_start(arglist, format);
  PCU_Debug_Printv_2(&global_pcu, format, arglist);
  va_end(arglist);
}

/** \brief Similar to PCU_Comm_Sender, returns the rank as an argument. */
int PCU_Comm_From(int* from_rank)
{
  return PCU_Comm_From_2(&global_pcu, from_rank);
}

/** \brief Returns in * \a size the bytes in the current received buffer
  \details This function should be called after a successful PCU_Comm_Receive.
  The size returned will be the total received size regardless of how
  much unpacking has been done.
 */
int PCU_Comm_Received(size_t* size)
{
  return PCU_Comm_Received_2(&global_pcu, size);
}

/** \brief Extracts a block of data from the current received buffer.
  \details This function should be called after a successful PCU_Comm_Receive.
  The next \a size bytes of the current received buffer are unpacked,
  and an internal pointer to that data is returned.
  The returned pointer must not be freed by the user.
 */
void* PCU_Comm_Extract(size_t size)
{
  return PCU_Comm_Extract_2(&global_pcu, size);
}

/** \brief Reinitializes PCU with a new MPI communicator.
 \details All of PCU's logic is based off two duplicates
 of this communicator, so you can safely get PCU to act
 on sub-groups of processes using this function.
 This call should be collective over all processes
 in the previous communicator. This is a very heavy weight function
 and should be used sparingly.
 */
void PCU_Switch_Comm(MPI_Comm new_comm)
{
  PCU_Switch_Comm_2(&global_pcu, new_comm);
}

/** \brief Return the current MPI communicator
  \details Returns the communicator given to the
  most recent PCU_Switch_Comm call, or MPI_COMM_WORLD
  otherwise.
 */
MPI_Comm PCU_Get_Comm(void)
{
  PCU_Get_Comm_2(&global_pcu);
}

/** \brief Return the time in seconds since some time in the past
 */
double PCU_Time(void)
{
  return PCU_Time_2(&global_pcu);
}

void PCU_Protect(void)
{
  PCU_Protect_2(&global_pcu);
}
