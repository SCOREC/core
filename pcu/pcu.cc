#include "PCU.h"
#include "PCUObj.h"
#include "reel.h"
#include <cstdarg>



extern "C" {

int PCU_Comm_Init(PCU_t* h) {
  if (h->ptr != nullptr)
    reel_fail("nested calls to Comm_Init");
  pcu::PCU* pcu_object = new pcu::PCU(MPI_COMM_WORLD);
  h->ptr = static_cast<void*>(pcu_object);
  return PCU_SUCCESS;
}

int PCU_Comm_Free(PCU_t* h) {
  if (h->ptr == nullptr)
    reel_fail("Comm_Free called before Comm_Init");
  delete static_cast<pcu::PCU*>(h->ptr);
  h->ptr = nullptr;
  return PCU_SUCCESS;
}

int PCU_Comm_Self(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Self called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Self();
}

int PCU_Comm_Peers(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Peers called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Peers();
}

void PCU_Comm_Begin(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Begin called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Begin();
}

/** \brief Packs data to be sent to \a to_rank.
  \details This function appends the block of \a size bytes starting
  at \a data to the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Pack(PCU_t h, int to_rank, const void *data, size_t size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Pack called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Pack(to_rank, data, size);
}

/** \brief Sends all buffers for this communication phase.
  \details This function should be called by all threads in the MPI job
  after calls to PCU_Comm_Pack or PCU_Comm_Write and before calls
  to PCU_Comm_Listen or PCU_Comm_Read.
  All buffers from this thread are sent out and receiving
  may begin after this call.
 */
int PCU_Comm_Send(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Send called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Send();
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
bool PCU_Comm_Listen(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Listen called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Listen();
}

/** \brief Returns in * \a from_rank the sender of the current received buffer.
  \details This function should be called after a successful PCU_Comm_Listen.
 */
int PCU_Comm_Sender(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Sender called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Sender();
}

bool PCU_Comm_Unpacked(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Unpacked called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Unpacked();
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
int PCU_Comm_Unpack(PCU_t h, void *data, size_t size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Unpack called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Unpack(data, size);
}

void PCU_Comm_Order(PCU_t h, bool on) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Order called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Order(on);
}

/** \brief Blocking barrier over all threads. */
void PCU_Barrier(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Barrier called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Barrier();
}

/** \brief Performs an Allreduce sum of double arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n doubles.
  After this call, p[i] will contain the sum of all p[i]'s
  given by each rank.
  */
void PCU_Add_Doubles(PCU_t h, double *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Add_Doubles called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Add(p, n);
}

double PCU_Add_Double(PCU_t h, double x) {
  if (h.ptr == nullptr)
    reel_fail("Add_Double called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Add(x);
}

/** \brief Performs an Allreduce minimum of double arrays.
 */
void PCU_Min_Doubles(PCU_t h, double *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Min_Doubles called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Min(p, n);
}

double PCU_Min_Double(PCU_t h, double x) {
  if (h.ptr == nullptr)
    reel_fail("Min_Double called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Min(x);
}

/** \brief Performs an Allreduce maximum of double arrays.
 */
void PCU_Max_Doubles(PCU_t h, double *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Max_Doubles called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Max(p, n);
}

double PCU_Max_Double(PCU_t h, double x) {
  if (h.ptr == nullptr)
    reel_fail("Max_Double called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Max(x);
}

/** \brief Performs an Allreduce sum of integers
 */
void PCU_Add_Ints(PCU_t h, int *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Add_Ints called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Add(p, n);
}

int PCU_Add_Int(PCU_t h, int x) {
  if (h.ptr == nullptr)
    reel_fail("Add_Int called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Add(x);
}

/** \brief Performs an Allreduce sum of long integers
 */
void PCU_Add_Longs(PCU_t h, long *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Add_Longs called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Add(p, n);
}

long PCU_Add_Long(PCU_t h, long x) {
  if (h.ptr == nullptr)
    reel_fail("Add_Long called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Add(x);
}

/** \brief Performs an Allreduce sum of size_t unsigned integers
 */
void PCU_Add_SizeTs(PCU_t h, size_t *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Add_SizeTs called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Add(p, n);
}

size_t PCU_Add_SizeT(PCU_t h, size_t x) {
  if (h.ptr == nullptr)
    reel_fail("Add_SizeT called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Add(x);
}

/** \brief Performs an Allreduce minimum of size_t unsigned integers
 */
void PCU_Min_SizeTs(PCU_t h, size_t *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Min_SizeTs called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Min(p, n);
}

size_t PCU_Min_SizeT(PCU_t h, size_t x) {
  if (h.ptr == nullptr)
    reel_fail("Min_SizeT called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Min(x);
}

/** \brief Performs an Allreduce maximum of size_t unsigned integers
 */
void PCU_Max_SizeTs(PCU_t h, size_t *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Max_SizeTs called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Max(p, n);
}

size_t PCU_Max_SizeT(PCU_t h, size_t x) {
  if (h.ptr == nullptr)
    reel_fail("Max_SizeT called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Max(x);
}

/** \brief Performs an exclusive prefix sum of integer arrays.
  \details This function must be called by all ranks at
  the same time. \a p must point to an array of \a n integers.
  After this call, p[i] will contain the sum of all p[i]'s
  given by ranks lower than the calling rank.
  */
void PCU_Exscan_Ints(PCU_t h, int *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Exscan_Ints called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Exscan(p, n);
}

int PCU_Exscan_Int(PCU_t h, int x) {
  if (h.ptr == nullptr)
    reel_fail("Exscan_Int called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Exscan(x);
}

/** \brief See PCU_Exscan_Ints */
void PCU_Exscan_Longs(PCU_t h, long *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Exscan_Longs called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Exscan(p, n);
}

long PCU_Exscan_Long(PCU_t h, long x) {
  if (h.ptr == nullptr)
    reel_fail("Exscan_Long called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Exscan(x);
}

/** \brief Performs an Allreduce minimum of int arrays.
 */
void PCU_Min_Ints(PCU_t h, int *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Min_Ints called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Min(p, n);
}

int PCU_Min_Int(PCU_t h, int x) {
  if (h.ptr == nullptr)
    reel_fail("Min_Int called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Min(x);
}

/** \brief Performs an Allreduce maximum of int arrays.
 */
void PCU_Max_Ints(PCU_t h, int *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Max_Ints called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Max(p, n);
}

int PCU_Max_Int(PCU_t h, int x) {
  if (h.ptr == nullptr)
    reel_fail("Max_Int called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Max(x);
}
/** \brief Performs an Allreduce maximum of long arrays.
 */
void PCU_Max_Longs(PCU_t h, long *p, size_t n) {
  if (h.ptr == nullptr)
    reel_fail("Max_Longs called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->Max(p, n);
}

long PCU_Max_Long(PCU_t h, long x) {
  if (h.ptr == nullptr)
    reel_fail("Max_Long called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Max(x);
}

/** \brief Performs a parallel logical OR reduction
 */
int PCU_Or(PCU_t h, int c) {
  if (h.ptr == nullptr)
    reel_fail("Or called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Or(c);
}

/** \brief Performs a parallel logical AND reduction
 */
int PCU_And(PCU_t h, int c) {
  if (h.ptr == nullptr)
    reel_fail("And called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->And(c);
}

/** \brief Returns the unique rank of the calling process.
 */
int PCU_Proc_Self(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Proc_Self called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Self();
}

/** \brief Returns the number of processes.
 */
int PCU_Proc_Peers(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Proc_Peers called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Peers();
}

/** \brief Similar to PCU_Comm_Self, returns the rank as an argument.
 */
int PCU_Comm_Rank(PCU_t h, int *rank) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Rank called before Comm_Init");
  *rank = static_cast<pcu::PCU*>(h.ptr)->Self();
  return PCU_SUCCESS;
}

/** \brief Similar to PCU_Comm_Peers, returns the size as an argument. */
int PCU_Comm_Size(PCU_t h, int *size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Size called before Comm_Init");
  *size = static_cast<pcu::PCU*>(h.ptr)->Peers();
  return PCU_SUCCESS;
}

/** \brief Returns true iff PCU has been initialized */
bool PCU_Comm_Initialized(PCU_t h) { return h.ptr != nullptr; }



/** \brief Returns in * \a size the number of bytes being sent to \a to_rank.
  \details Returns the size of the buffer being sent to \a to_rank.
  This function should be called after PCU_Comm_Start and before
  PCU_Comm_Send.
 */
int PCU_Comm_Packed(PCU_t h, int to_rank, size_t *size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Packed called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Packed(to_rank, size);
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
int PCU_Comm_Write(PCU_t h, int to_rank, const void *data, size_t size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Write called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Write(to_rank, data, size);
}

/** \brief Convenience wrapper over Listen and Unpacked */
bool PCU_Comm_Receive(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Receive called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Receive();
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
bool PCU_Comm_Read(PCU_t h, int *from_rank, void **data, size_t *size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Read called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Read(from_rank, data, size);
}

void PCU_Debug_Open(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Debug_Open called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->DebugOpen();
}

/** \brief like fprintf, contents go to debugN.txt */
void PCU_Debug_Print(PCU_t h, const char *format, ...) {
  if (h.ptr == nullptr)
    reel_fail("Debug_Print called before Comm_Init");
  va_list arglist;
  va_start(arglist, format);
  static_cast<pcu::PCU*>(h.ptr)->DebugPrint(format, arglist);
  va_end(arglist);
}

/** \brief Similar to PCU_Comm_Sender, returns the rank as an argument. */
int PCU_Comm_From(PCU_t h, int *from_rank) {
  if (h.ptr == nullptr)
    reel_fail("Comm_From called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->From(from_rank);
}

/** \brief Returns in * \a size the bytes in the current received buffer
  \details This function should be called after a successful PCU_Comm_Receive.
  The size returned will be the total received size regardless of how
  much unpacking has been done.
 */
int PCU_Comm_Received(PCU_t h, size_t *size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Received called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Received(size);
}

/** \brief Extracts a block of data from the current received buffer.
  \details This function should be called after a successful PCU_Comm_Receive.
  The next \a size bytes of the current received buffer are unpacked,
  and an internal pointer to that data is returned.
  The returned pointer must not be freed by the user.
 */
void *PCU_Comm_Extract(PCU_t h, size_t size) {
  if (h.ptr == nullptr)
    reel_fail("Comm_Extract called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->Extract(size);
}

/** \brief Reinitializes PCU with a new MPI communicator.
 \details All of PCU's logic is based off two duplicates
 of this communicator, so you can safely get PCU to act
 on sub-groups of processes using this function.
 This call should be collective over all processes
 in the previous communicator. This is a very heavy weight function
 and should be used sparingly.
 */
void PCU_Switch_Comm(PCU_t h, MPI_Comm new_comm) {
  if (h.ptr == nullptr)
    reel_fail("Switch_Comm called before Comm_Init");
  static_cast<pcu::PCU*>(h.ptr)->SwitchMPIComm(new_comm);
}

/** \brief Return the current MPI communicator
  \details Returns the communicator given to the
  most recent PCU_Switch_Comm call, or MPI_COMM_WORLD
  otherwise.
 */
MPI_Comm PCU_Get_Comm(PCU_t h) {
  if (h.ptr == nullptr)
    reel_fail("Get_Comm called before Comm_Init");
  return static_cast<pcu::PCU*>(h.ptr)->GetMPIComm();
}

/** \brief Return the time in seconds since some time in the past
 */
double PCU_Time(void) { return pcu::Time(); }

void PCU_Protect(void) { return pcu::Protect(); }

double PCU_GetMem(void) { return pcu::GetMem(); }


}