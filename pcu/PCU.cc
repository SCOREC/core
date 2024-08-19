#include "PCU.h"
#include "noto_malloc.h"
#include "pcu_mem.h"
#include "pcu_mpi.h"
#include "pcu_msg.h"
#include "pcu_order.h"
#include "reel.h"
#include <algorithm>
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <climits>
#include <cstring>
#include <cerrno>
#include <cstdarg>
namespace pcu {

int PCU::Peers() const noexcept { return pcu_mpi_size(mpi_); }
int PCU::Self() const noexcept { return pcu_mpi_rank(mpi_); }
void PCU::Begin() noexcept { pcu_msg_start(mpi_, msg_); }
int PCU::Pack(int to_rank, const void *data, size_t size) noexcept {
  if ((to_rank < 0) || (to_rank >= Peers()))
    reel_fail("Invalid rank in Comm_Pack");
  if (size > (size_t)INT_MAX) {
    fprintf(stderr, "ERROR Attempting to pack a PCU message whose size exceeds "
                    "INT_MAX... exiting\n");
    abort();
  }
  memcpy(pcu_msg_pack(msg_, to_rank, size), data, size);
  return PCU_SUCCESS;
}

int PCU::Send() noexcept {
  pcu_msg_send(mpi_, msg_);
  return PCU_SUCCESS;
}
bool PCU::Receive() noexcept {
  while (Unpacked())
    if (!Listen())
      return false;
  return true;
}
bool PCU::Listen() noexcept {
  if (msg_->order)
    return pcu_order_receive(mpi_, msg_->order, msg_);
  return pcu_msg_receive(mpi_, msg_);
}
int PCU::Sender() noexcept {
  if (msg_->order)
    return pcu_order_received_from(msg_->order);
  return pcu_msg_received_from(msg_);
}
bool PCU::Unpacked() noexcept {
  if (msg_->order)
    return pcu_order_unpacked(msg_->order);
  return pcu_msg_unpacked(msg_);
}
int PCU::Unpack(void *data, size_t size) noexcept {
  if (msg_->order)
    memcpy(data, pcu_order_unpack(msg_->order, size), size);
  else
    memcpy(data, pcu_msg_unpack(msg_, size), size);
  return PCU_SUCCESS;
}

int PCU::Write(int to_rank, const void *data, size_t size) noexcept {
  if ((to_rank < 0) || (to_rank >= Peers()))
    reel_fail("Invalid rank in Comm_Write");
  PCU_MSG_PACK(msg_, to_rank, size);
  memcpy(pcu_msg_pack(msg_, to_rank, size), data, size);
  return PCU_SUCCESS;
}
bool PCU::Read(int *from_rank, void **data, size_t *size) noexcept {
  if (!Receive())
    return false;
  *from_rank = Sender();
  Unpack(size, sizeof(*size));
  *data = Extract(*size);
  return true;
}
void PCU::Order(bool on) {
  if (on && (!msg_->order))
    msg_->order = pcu_order_new();
  if ((!on) && msg_->order) {
    pcu_order_free(msg_->order);
    msg_->order = NULL;
  }
}
void PCU::Barrier() { pcu_barrier(mpi_, &(msg_->coll)); }
int PCU::Or(int c) noexcept { return Max(c); }
int PCU::And(int c) noexcept { return Min(c); }
int PCU::Packed(int to_rank, size_t *size) noexcept {
  if ((to_rank < 0) || (to_rank >= Peers()))
    reel_fail("Invalid rank in Comm_Packed");
  *size = pcu_msg_packed(msg_, to_rank);
  return PCU_SUCCESS;
}
int PCU::From(int *from_rank) noexcept {
  if (msg_->order)
    *from_rank = pcu_order_received_from(msg_->order);
  else
    *from_rank = pcu_msg_received_from(msg_);
  return PCU_SUCCESS;
}
int PCU::Received(size_t *size) noexcept {
  if (msg_->order)
    *size = pcu_order_received_size(msg_->order);
  else
    *size = pcu_msg_received_size(msg_);
  return PCU_SUCCESS;
}
void *PCU::Extract(size_t size) noexcept {
  if (msg_->order)
    return pcu_order_unpack(msg_->order, size);
  return pcu_msg_unpack(msg_, size);
}

static void safe_mkdir(const char *path, mode_t mode) {
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
    reel_fail("PCU: could not create directory \"%s\"\n", path);
}

static void append(char *s, size_t size, const char *format, ...) {
  int len = strlen(s);
  va_list ap;
  va_start(ap, format);
  vsnprintf(s + len, size - len, format, ap);
  va_end(ap);
}

void PCU::DebugOpen() noexcept {
  const int fanout = 2048;
  const int bufsize = 1024;
  char *path = (char *)noto_malloc(bufsize);
  path[0] = '\0';
  if (Peers() > fanout) {
    mode_t const dir_perm = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
    strcpy(path, "debug/");
    safe_mkdir(path, dir_perm);
    int self = Self();
    append(path, bufsize, "%d/", self / fanout);
    if (self % fanout == 0)
      safe_mkdir(path, dir_perm);
    Barrier();
  }

  append(path, bufsize, "%s", "debug");
  if (!msg_->file)
    msg_->file = pcu_open_parallel(GetCHandle(), path, "txt");
  noto_free(path);
}



double GetMem() noexcept { return pcu_get_mem(); }
void Protect() noexcept { reel_protect(); }
double Time() noexcept { return MPI_Wtime(); }

void PCU::DebugPrint(const char *format, ...) noexcept {
  va_list args;
  va_start(args, format);
  DebugPrint(format, args);
  va_end(args);
}
void PCU::DebugPrint(const char *format, va_list args) noexcept {
  if (!msg_->file)
    return; // Print is a no-op if no file is open
  vfprintf(msg_->file, format, args);
  fflush(msg_->file);
}
PCU::PCU(MPI_Comm comm) {
  mpi_ = new pcu_mpi_t;
  msg_ = new pcu_msg;
  pcu_mpi_init(comm, mpi_);
  pcu_make_msg(msg_);
  /* turn ordering on by default, call
     PCU_Comm_Order(false) after PCU_Comm_Init
     to disable this */
  Order(true);
}
PCU::~PCU() noexcept {
  pcu_mpi_finalize(mpi_);
  delete mpi_;
  pcu_free_msg(msg_);
  delete msg_;
}
PCU::PCU(PCU &&other) noexcept {
  std::swap(mpi_, other.mpi_);
  std::swap(msg_, other.msg_);
}
PCU &PCU::operator=(PCU && other) noexcept {
  std::swap(mpi_, other.mpi_);
  std::swap(msg_, other.msg_);
  return *this;
}
MPI_Comm PCU::GetMPIComm() const noexcept { return mpi_->original_comm; }

MPI_Comm PCU::SwitchMPIComm(MPI_Comm newcomm) noexcept {
  if(newcomm == mpi_->original_comm) {
    return mpi_->original_comm;
  }
  auto original_comm = mpi_->original_comm;
  pcu_mpi_finalize(mpi_);
  pcu_mpi_init(newcomm, mpi_);
  return original_comm;
}

/* template implementations */
template <typename T> void PCU::Add(T *p, size_t n) noexcept {
  pcu_allreduce(
      mpi_, &(msg_->coll),
      [](void *local, void *incoming, size_t size) {
        auto *a = static_cast<T *>(local);
        auto *b = static_cast<T *>(incoming);
        size_t n = size / sizeof(T);
        for (size_t i = 0; i < n; ++i)
          a[i] += b[i];
      },
      p, n * sizeof(T));
}
template <typename T> T PCU::Add(T p) noexcept {
  Add(&p, 1);
  return p;
}
template <typename T> void PCU::Min(T *p, size_t n) noexcept {
  pcu_allreduce(
      mpi_, &(msg_->coll),
      [](void *local, void *incoming, size_t size) {
        auto *a = static_cast<T *>(local);
        auto *b = static_cast<T *>(incoming);
        size_t n = size / sizeof(T);
        for (size_t i = 0; i < n; ++i)
          a[i] = std::min(a[i], b[i]);
      },
      p, n * sizeof(T));
}
template <typename T> T PCU::Min(T p) noexcept {
  Min(&p, 1);
  return p;
}
template <typename T> void PCU::Max(T *p, size_t n) noexcept {
  pcu_allreduce(
      mpi_, &(msg_->coll),
      [](void *local, void *incoming, size_t size) {
        auto *a = static_cast<T *>(local);
        auto *b = static_cast<T *>(incoming);
        size_t n = size / sizeof(T);
        for (size_t i = 0; i < n; ++i)
          a[i] = std::max(a[i], b[i]);
      },
      p, n * sizeof(T));
}
template <typename T> T PCU::Max(T p) noexcept {
  Max(&p, 1);
  return p;
}
template <typename T> void PCU::Exscan(T *p, size_t n) noexcept {
  auto *originals = (T *)noto_malloc(sizeof(T) * n);
  for (size_t i = 0; i < n; ++i)
    originals[i] = p[i];
  pcu_scan(
      mpi_, &(msg_->coll),
      [](void *local, void *incoming, size_t size) {
        auto *a = static_cast<T *>(local);
        auto *b = static_cast<T *>(incoming);
        size_t n = size / sizeof(T);
        for (size_t i = 0; i < n; ++i)
          a[i] += b[i];
      },
      p, n * sizeof(T));
  // convert inclusive scan to exclusive
  for (size_t i = 0; i < n; ++i)
    p[i] -= originals[i];
  noto_free(originals);
}
template <typename T> T PCU::Exscan(T p) noexcept {
  Exscan(&p, 1);
  return p;
}
#define PCU_EXPL_INST_DECL(T)                                                  \
  template void PCU::Add<T>(T * p, size_t n) noexcept;                         \
  template T PCU::Add<T>(T p) noexcept;                                        \
  template void PCU::Min<T>(T * p, size_t n) noexcept;                         \
  template T PCU::Min<T>(T p) noexcept;                                        \
  template void PCU::Max<T>(T * p, size_t n) noexcept;                         \
  template T PCU::Max<T>(T p) noexcept;                                        \
  template void PCU::Exscan<T>(T * p, size_t n) noexcept;                      \
  template T PCU::Exscan<T>(T p) noexcept;
PCU_EXPL_INST_DECL(int)
PCU_EXPL_INST_DECL(size_t)
PCU_EXPL_INST_DECL(long)
PCU_EXPL_INST_DECL(double)
#undef PCU_EXPL_INST_DECL

} // namespace pcu
