/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_H
#define PCU_H

#define PCU_SUCCESS 0
#define PCU_FAILURE -1

#include <mpi.h>

#ifdef __cplusplus
#include <cstddef>
#include <cstdio>
extern "C" {
#else
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#endif

/*library init/finalize*/
int PCU_Comm_Init(void);
int PCU_Comm_Free(void);

/*rank/size functions*/
int PCU_Comm_Self(void);
int PCU_Comm_Peers(void);

/*recommended message passing API*/
void PCU_Comm_Begin(void);
int PCU_Comm_Pack(int to_rank, const void* data, size_t size);
#define PCU_COMM_PACK(to_rank,object)\
PCU_Comm_Pack(to_rank,&(object),sizeof(object))
int PCU_Comm_Send(void);
bool PCU_Comm_Receive(void);
bool PCU_Comm_Listen(void);
int PCU_Comm_Sender(void);
bool PCU_Comm_Unpacked(void);
int PCU_Comm_Unpack(void* data, size_t size);
#define PCU_COMM_UNPACK(object)\
PCU_Comm_Unpack(&(object),sizeof(object))

/*turns deterministic ordering for the
  above API on/off*/
void PCU_Comm_Order(bool on);

/*collective operations*/
void PCU_Barrier(void);
void PCU_Add_Doubles(double* p, size_t n);
void PCU_Min_Doubles(double* p, size_t n);
void PCU_Max_Doubles(double* p, size_t n);
void PCU_Add_Ints(int* p, size_t n);
void PCU_Add_Longs(long* p, size_t n);
void PCU_Exscan_Ints(int* p, size_t n);
void PCU_Exscan_Longs(long* p, size_t n);
void PCU_Min_Ints(int* p, size_t n);
void PCU_Max_Ints(int* p, size_t n);
int PCU_Or(int c);

/*thread functions*/
typedef void* (*PCU_Thrd_Func)(void*);
int PCU_Thrd_Run(int nthreads, PCU_Thrd_Func function, void** in_out);
int PCU_Thrd_Self(void);
int PCU_Thrd_Peers(void);
void PCU_Thrd_Barrier(void);
void PCU_Thrd_Lock(void);
void PCU_Thrd_Unlock(void);

/*process-level self/peers (mpi wrappers)*/
int PCU_Proc_Self(void);
int PCU_Proc_Peers(void);

/*IPComMan replacement API*/
int PCU_Comm_Write(int to_rank, const void* data, size_t size);
#define PCU_COMM_WRITE(to,data) \
PCU_Comm_Write(to,&(data),sizeof(data))
bool PCU_Comm_Read(int* from_rank, void** data, size_t* size);

/*Debug file I/O API*/
void PCU_Debug_Open(void);
#ifdef __GNUC__
void PCU_Debug_Print(const char* format, ...)
  __attribute__((format(printf,1,2)));
#else
void PCU_Debug_Print(const char* format, ...);
#endif

/*lesser-used APIs*/
bool PCU_Comm_Initialized(void);
int PCU_Comm_Packed(int to_rank, size_t* size);
int PCU_Comm_From(int* from_rank);
int PCU_Comm_Received(size_t* size);
void* PCU_Comm_Extract(size_t size);
int PCU_Comm_Rank(int* rank);
int PCU_Comm_Size(int* size);

/*deprecated method enum*/
#ifdef __cplusplus
enum PCU_Method { PCU_GLOBAL_METHOD, PCU_LOCAL_METHOD };
#else
typedef enum { PCU_GLOBAL_METHOD, PCU_LOCAL_METHOD } PCU_Method;
#endif
int PCU_Comm_Start(PCU_Method method);

/*special MPI_Comm replacement API*/
void PCU_Switch_Comm(MPI_Comm new_comm);
MPI_Comm PCU_Get_Comm(void);

/*stack trace helpers using GNU/Linux*/
void PCU_Trace(void);
void PCU_Protect(void);

/*MPI_Wtime() equivalent*/
double PCU_Time(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
