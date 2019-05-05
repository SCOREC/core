%module pcu
%{
#include <mpi.h>
#include <cstddef>
#include <cstdio>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <PCU.h>
#include <pcu_util.h>
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);


/* ==== FROM PCU.h ====*/
#define PCU_SUCCESS 0
#define PCU_FAILURE -1

/*library init/finalize*/
int PCU_Comm_Init(void);
int PCU_Comm_Free(void);

int PCU_Comm_Self(void);
int PCU_Comm_Peers(void)
;
double PCU_Time(void); 

/* ==== FROM pcu_util.h ====*/
void PCU_Assert_Fail(const char* msg);

/* This are defined as macros in the .h file. Apparaently, it is OK to Lie
to SWIG that these are functions ;) 
*/
void PCU_ALWAYS_ASSERT(int cond);
void PCU_ALWAYS_ASSERT_VERBOSE(int cond, const char* msg);
