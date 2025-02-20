#ifndef SCOREC_PCU_H
#define SCOREC_PCU_H

#include <cstdlib>
#include <cstdarg> //va_list
#include "pcu_defines.h"

struct pcu_msg_struct;
struct pcu_mpi_struct;

namespace pcu {
class PCU {
public:
  PCU();
  explicit PCU(PCU_Comm comm);
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
  [[deprecated("Use PCU::GetComm instead.")]]
  [[nodiscard]] PCU_Comm GetMPIComm() const noexcept;
  /**
   * @brief Get the underlying communicator which may be an MPI_Comm.
   */
  [[nodiscard]] PCU_Comm GetComm() const noexcept;
  /** @brief Check if the original PCU_Comm is owned by this object.
   * If true, it will be freed during destruction. */
  bool OwnsComm() const noexcept;
  /** @brief Set ownership of the orignal PCU_Comm passed to the constructor.
   *
   * This function can enable or disable ownership. If a communicator created
   * with PCU::Split() is disowned, it should be freed by the user.
   *
   * @param on true to enable ownership.
   */
  void OwnsComm(bool on) const noexcept;

  [[nodiscard]] PCU_t GetCHandle() {PCU_t h; h.ptr=this; return h;}
  /*recommended message passing API*/
  void Begin() noexcept;
  int Pack(int to_rank, const void *data, size_t size) noexcept;
  template<typename T> int Pack(int to_rank, T& data) noexcept {
    return Pack(to_rank, &(data), sizeof(data));
  }
  template<typename T> int Pack(int to_rank, T*& data) noexcept {
    return Pack(to_rank, &(data), sizeof(data));
  }

  int Send() noexcept;
  bool Receive() noexcept;
  bool Listen() noexcept;
  int Sender() noexcept;
  bool Unpacked() noexcept;
  int Unpack(void *data, size_t size) noexcept;
  template<typename T> int Unpack(T& data) noexcept {
    return Unpack(&(data), sizeof(data));
  }
  template<typename T> int Unpack(T*& data) noexcept {
    return Unpack(&(data), sizeof(data));
  }
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
  template <typename T> void Allgather(T *send, T *recv, size_t n) noexcept;

  /*bitwise operations*/
  [[nodiscard]] int Or(int c) noexcept;
  [[nodiscard]] int And(int c) noexcept;

  /**
   * @brief Split a communicator into distinct subgroups.
   *
   * The resulting communicator is marked owned and automatically free the
   * underlying communicator. This can be disabled with PCU::OwnsComm(bool). In
   * that case, the user is responsible for cleanup.
   *
   * @param color subgroup indicator.
   * @param key used for subgroup ordering; specify 0 if you don't care.
   * @return a new communicator defined on the resulting subgroup.
   */
  PCU* Split(int color, int key) noexcept;

  /*lesser-used APIs*/
  int Packed(int to_rank, size_t *size) noexcept;
  int From(int *from_rank) noexcept;
  int Received(size_t *size) noexcept;
  void *Extract(size_t size) noexcept;

  void DebugPrint(const char* format, ...) noexcept PCU_FORMAT_ATTRIBUTE(2, 3)
  void DebugPrint(const char* format, va_list args) noexcept;
  /* Debug functions */
  void DebugOpen() noexcept;

  [[deprecated("Use PCU::SwitchComm instead.")]]
#if __cplusplus >= 201703L
  [[nodiscard]]
#endif
  PCU_Comm SwitchMPIComm(PCU_Comm) noexcept;

#if __cplusplus >= 201703L
  [[nodiscard]]
#endif
  PCU_Comm SwitchComm(PCU_Comm) noexcept;

private:
  pcu_msg_struct *msg_;
  pcu_mpi_struct *mpi_;
};
/*stack trace helpers using GNU/Linux*/
void Protect() noexcept;
/*Memory usage*/
[[nodiscard]] double GetMem() noexcept;
/*MPI_Wtime() equivalent*/
[[nodiscard]] double Time() noexcept;

PCU* PCU_GetGlobal();

/* explicit instantiations of template functions */
#define PCU_EXPL_INST_DECL(T)                                                  \
  extern template void PCU::Add<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Add<T>(T p) noexcept;                                 \
  extern template void PCU::Min<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Min<T>(T p) noexcept;                                 \
  extern template void PCU::Max<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Max<T>(T p) noexcept;                                 \
  extern template void PCU::Exscan<T>(T * p, size_t n) noexcept;               \
  extern template T PCU::Exscan<T>(T p) noexcept;                              \
  extern template void PCU::Allgather<T>(T *send, T *recv, size_t n) noexcept;
PCU_EXPL_INST_DECL(int)
PCU_EXPL_INST_DECL(size_t)
PCU_EXPL_INST_DECL(long)
PCU_EXPL_INST_DECL(double)
#undef PCU_EXPL_INST_DECL

} // namespace pcu

#endif // PCUOBJ_H


