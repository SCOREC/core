#ifndef SCOREC_PCU_PCU_DEFINES_H
#define SCOREC_PCU_PCU_DEFINES_H

#include <SCOREC_config.h>
#ifndef PUMI_NO_MPI
#include <mpi.h>
#endif

#define PCU_SUCCESS 0
#define PCU_FAILURE -1
#ifdef __GNUC__
#define PCU_FORMAT_ATTRIBUTE(...) \
  __attribute__((format(printf, ##__VA_ARGS__)));
#else
#define PCU_FORMAT_ATTRIBUTE(format, ...)
#endif

#ifdef __cplusplus
extern "C"{
#endif

#ifndef PUMI_NO_MPI
typedef MPI_Comm PCU_Comm;
typedef MPI_Request PCU_Request;
#define PCU_ANY_SOURCE MPI_ANY_SOURCE
#else
typedef int PCU_Comm;
typedef int PCU_Request;
#define PCU_ANY_SOURCE -1
#endif

struct PCU_t {
    void* ptr;
};

#ifdef __cplusplus
}
#endif

#endif // SCOREC_PCU_PCU_DEFINES_H
