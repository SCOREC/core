#ifndef SCOREC_PCU_H
#define SCOREC_PCU_H

#include <memory>
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
  template <typename T>
  void Allgather(const T *send, T *recv, size_t n) noexcept;

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
  std::unique_ptr<PCU> Split(int color, int key) noexcept;

  /**
   * @brief Duplicate the underlying communicator.
   *
   * It is the user's responsiblity to cleanup the returned communicator. This
   * function may be used to initialize other libraries using the PCU-defined
   * communication group.
   *
   * If SCOREC::core was compiled with the SCOREC_NO_MPI flag, the return value
   * is not meaningful.
   *
   * @param[out] newcomm The output address for the new communicator copy.
   * @return 0 on success.
   */
  int DupComm(PCU_Comm* newcomm) const noexcept;

  /*lesser-used APIs*/
  int Packed(int to_rank, size_t *size) noexcept;
  int From(int *from_rank) noexcept;
  int Received(size_t *size) noexcept;
  void *Extract(size_t size) noexcept;

  /*
   * Debug print function with printf-style formatting
   * This function is excluded from SWIG parsing due to complex macro usage
   * with variadic argument. Also, such functionalities is not needed for the Python API.
   */
  #ifndef SWIG
  void DebugPrint(const char* format, ...) noexcept PCU_FORMAT_ATTRIBUTE(2, 3)
  void DebugPrint(const char* format, va_list args) noexcept;
  #endif // SWIG
  /* Debug functions */
  void DebugOpen() noexcept;

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

/**
 * @brief Initialize the underlying parallel library.
 *
 * This may be MPI (or a stub, given SCOREC_NO_MPI). This function abstracts
 * the difference.
 */
void Init(int *argc, char ***argv);
/**
 * @brief Finalize the underlying parallel library.
 *
 * This may be MPI (or a stub, given SCOREC_NO_MPI). This function abstracts
 * the difference.
 */
void Finalize();

/*
 * Explicit instantiations of template functions,
 * ignored by SWIG to avoid difficulties when parsing macros involving
 * template functions. Templates are initialized manually in the interface file.
 */
#ifndef SWIG
#define PCU_EXPL_INST_DECL(T)                                                  \
  extern template void PCU::Add<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Add<T>(T p) noexcept;                                 \
  extern template void PCU::Min<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Min<T>(T p) noexcept;                                 \
  extern template void PCU::Max<T>(T * p, size_t n) noexcept;                  \
  extern template T PCU::Max<T>(T p) noexcept;                                 \
  extern template void PCU::Exscan<T>(T * p, size_t n) noexcept;               \
  extern template T PCU::Exscan<T>(T p) noexcept;                              \
  extern template                                                              \
  void PCU::Allgather<T>(const T *send, T *recv, size_t n) noexcept;
PCU_EXPL_INST_DECL(int)
PCU_EXPL_INST_DECL(size_t)
PCU_EXPL_INST_DECL(long)
PCU_EXPL_INST_DECL(double)
#undef PCU_EXPL_INST_DECL
#endif // SWIG

} // namespace pcu

#endif // PCUOBJ_H


