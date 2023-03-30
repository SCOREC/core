#ifndef SCOREC_PCU_PCUOBJ_H
#define SCOREC_PCU_PCUOBJ_H
#include <cstdlib>
#include <mpi.h>
#include "pcu_defines.h"

struct pcu_msg_struct;
struct pcu_mpi_struct;

namespace pcu {
class PCU {
public:
  explicit PCU(MPI_Comm comm);
  ~PCU() noexcept;
  PCU(PCU const &) = delete;
  PCU(PCU &&) noexcept;
  PCU &operator=(PCU const &) = delete;
  PCU &operator=(PCU &&) noexcept;
  /** @brief Returns the rank of the current process.
   *  @return The rank of the current process.
   */
  [[nodiscard]] int Self() const noexcept;
  /** @brief Returns the number of ranks in the communicator.
   *  @return The number of ranks in the communicator.
   */
  [[nodiscard]] int Peers() const noexcept;
  [[nodiscard]] MPI_Comm GetMPIComm() const noexcept;

  /*recommended message passing API*/
  void Begin() noexcept;
  int Pack(int to_rank, const void *data, size_t size) noexcept;
  int Send() noexcept;
  bool Receive() noexcept;
  bool Listen() noexcept;
  int Sender() noexcept;
  bool Unpacked() noexcept;
  int Unpack(void *data, size_t size) noexcept;
  /*IPComMan replacement API*/
  int Write(int to_rank, const void *data, size_t size) noexcept;
  bool Read(int *from_rank, void **data, size_t *size) noexcept;

  /*turns deterministic ordering for the
    above API on/off*/
  void Order(bool on);

  /*collective operations*/
  void Barrier();
  template <typename T> void Add(T *p, size_t n) noexcept;
  template <typename T> [[nodiscard]] T Add(T p) noexcept;
  template <typename T> void Min(T *p, size_t n) noexcept;
  template <typename T> [[nodiscard]] T Min(T p) noexcept;
  template <typename T> void Max(T *p, size_t n) noexcept;
  template <typename T> [[nodiscard]] T Max(T p) noexcept;
  template <typename T> void Exscan(T *p, size_t n) noexcept;
  template <typename T> [[nodiscard]] T Exscan(T p) noexcept;

  /*bitwise operations*/
  [[nodiscard]] int Or(int c) noexcept;
  [[nodiscard]] int And(int c) noexcept;

  /*lesser-used APIs*/
  int Packed(int to_rank, size_t *size) noexcept;
  int From(int *from_rank) noexcept;
  int Received(size_t *size) noexcept;
  void *Extract(size_t size) noexcept;

  void DebugPrint(const char* format, ...) noexcept PCU_FORMAT_ATTRIBUTE(2, 3);
  void DebugPrint(const char* format, va_list args) noexcept;
  /* Debug functions */
  void DebugOpen() noexcept;

private:
  pcu_msg_struct *msg_;
  pcu_mpi_struct *mpi_;
};
/*stack trace helpers using GNU/Linux*/
void Protect(void) noexcept;
/*Memory usage*/
[[nodiscard]] double GetMem(void) noexcept;
/*MPI_Wtime() equivalent*/
[[nodiscard]] double Time(void) noexcept;

/* explicit instantiations of template functions */
#define PCU_EXPL_INST_DECL(T)                                                  \
  extern template void PCU::Add<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Add<T>(T p) noexcept;                                 \
  extern template void PCU::Min<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Min<T>(T p) noexcept;                                 \
  extern template void PCU::Max<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Max<T>(T p) noexcept;                                 \
  extern template void PCU::Exscan<T>(T * p, size_t n) noexcept;               \
  extern template T PCU::Exscan<T>(T p) noexcept;
PCU_EXPL_INST_DECL(int)
PCU_EXPL_INST_DECL(size_t)
PCU_EXPL_INST_DECL(long)
PCU_EXPL_INST_DECL(double)
#undef PCU_EXPL_INST_DECL

} // namespace pcu
#undef PCU_FORMAT_ATTRIBUTE
#endif // SCOREC_PCU_PCUOBJ_H
