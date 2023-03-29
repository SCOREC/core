/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU2_H
#define PCU2_H

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

// FOR pcu_state_type ... use opaque ptr instead?
#include "pcu_msg.h"

typedef enum { pcu_state_uninit, pcu_state_init } pcu_state_enum ;
typedef struct {
  pcu_state_enum state;
  pcu_msg pmsg;
  pcu_mpi_t mpi;
} pcu_t;

/*library init/finalize*/
int PCU_Comm_Init_2(pcu_t* pcu);
int PCU_Comm_Free_2(pcu_t* pcu);

/*rank/size functions*/
int PCU_Comm_Self_2(pcu_t* pcu);
int PCU_Comm_Peers_2(pcu_t* pcu);

/*recommended message passing API*/
void PCU_Comm_Begin_2(pcu_t* pcu);
int PCU_Comm_Pack_2(pcu_t* pcu, int to_rank, const void* data, size_t size);
int PCU_Comm_Send_2(pcu_t* pcu);
bool PCU_Comm_Receive_2(pcu_t* pcu);
bool PCU_Comm_Listen_2(pcu_t* pcu);
int PCU_Comm_Sender_2(pcu_t* pcu);
bool PCU_Comm_Unpacked_2(pcu_t* pcu);
int PCU_Comm_Unpack_2(pcu_t* pcu, void* data, size_t size);

/*turns deterministic ordering for the
  above API on/off*/
void PCU_Comm_Order_2(pcu_t* pcu, bool on);

/*collective operations*/
void PCU_Barrier_2(pcu_t* pcu);
void PCU_Add_Doubles_2(pcu_t* pcu, double* p, size_t n);
double PCU_Add_Double_2(pcu_t* pcu, double x);
void PCU_Min_Doubles_2(pcu_t* pcu, double* p, size_t n);
double PCU_Min_Double_2(pcu_t* pcu, double x);
void PCU_Max_Doubles_2(pcu_t* pcu, double* p, size_t n);
double PCU_Max_Double_2(pcu_t* pcu, double x);
void PCU_Add_Ints_2(pcu_t* pcu, int* p, size_t n);
int PCU_Add_Int_2(pcu_t* pcu, int x);
void PCU_Add_Longs_2(pcu_t* pcu, long* p, size_t n);
long PCU_Add_Long_2(pcu_t* pcu, long x);
void PCU_Exscan_Ints_2(pcu_t* pcu, int* p, size_t n);
int PCU_Exscan_Int_2(pcu_t* pcu, int x);
void PCU_Exscan_Longs_2(pcu_t* pcu, long* p, size_t n);
long PCU_Exscan_Long_2(pcu_t* pcu, long x);
void PCU_Add_SizeTs_2(pcu_t* pcu, size_t* p, size_t n);
size_t PCU_Add_SizeT_2(pcu_t* pcu, size_t x);
void PCU_Min_SizeTs_2(pcu_t* pcu, size_t* p, size_t n);
size_t PCU_Min_SizeT_2(pcu_t* pcu, size_t x);
void PCU_Max_SizeTs_2(pcu_t* pcu, size_t* p, size_t n);
size_t PCU_Max_SizeT_2(pcu_t* pcu, size_t x);
void PCU_Min_Ints_2(pcu_t* pcu, int* p, size_t n);
int PCU_Min_Int_2(pcu_t* pcu, int x);
void PCU_Max_Ints_2(pcu_t* pcu, int* p, size_t n);
int PCU_Max_Int_2(pcu_t* pcu, int x);
void PCU_Max_Longs_2(pcu_t* pcu, long* p, size_t n);
long PCU_Max_Long_2(pcu_t* pcu, long x);
int PCU_Or_2(pcu_t* pcu, int c);
int PCU_And_2(pcu_t* pcu, int c);

/*process-level self/peers (mpi wrappers)*/
int PCU_Proc_Self_2(pcu_t* pcu);
int PCU_Proc_Peers_2(pcu_t* pcu);

/*IPComMan replacement API*/
int PCU_Comm_Write_2(pcu_t* pcu, int to_rank, const void* data, size_t size);
bool PCU_Comm_Read_2(pcu_t* pcu, int* from_rank, void** data, size_t* size);

/*Debug file I/O API*/
void PCU_Debug_Open_2(pcu_t* pcu);
void PCU_Debug_Printv_2(pcu_t* pcu, const char* format, va_list args);
#ifdef __GNUC__
void PCU_Debug_Print_2(pcu_t* pcu, const char* format, ...)
  __attribute__((format(printf,2,3)));
#else
void PCU_Debug_Print(pcu_t* pcu, const char* format, ...);
#endif

/*lesser-used APIs*/
bool PCU_Comm_Initialized_2(pcu_t* pcu);
int PCU_Comm_Packed_2(pcu_t* pcu, int to_rank, size_t* size);
int PCU_Comm_From_2(pcu_t* pcu, int* from_rank);
int PCU_Comm_Received_2(pcu_t* pcu, size_t* size);
void* PCU_Comm_Extract_2(pcu_t* pcu, size_t size);
int PCU_Comm_Rank_2(pcu_t* pcu, int* rank);
int PCU_Comm_Size_2(pcu_t* pcu, int* size);

int PCU_Comm_Start_2(pcu_t* pcu);

/*special MPI_Comm replacement API*/
void PCU_Switch_Comm_2(pcu_t* pcu, MPI_Comm new_comm);
MPI_Comm PCU_Get_Comm_2(pcu_t* pcu);

/*stack trace helpers using GNU/Linux*/
void PCU_Protect_2(pcu_t* pcu);

/*MPI_Wtime_() equivalent*/
double PCU_Time_2(pcu_t* pcu);

/*Memory usage*/
double PCU_GetMem_2(pcu_t* pcu);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
