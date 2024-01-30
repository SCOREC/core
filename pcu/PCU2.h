#ifndef PCU2_H
#define PCU2_H
#include "pcu_defines.h"
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

struct PCUHandle {
    void* ptr;
};

int PCU_Comm_Init2(PCUHandle* h);
int PCU_Comm_Free2(PCUHandle* h);

/*rank/size functions*/
int PCU_Comm_Self2(PCUHandle h);
int PCU_Comm_Peers2(PCUHandle h);

/*recommended message passing API*/
void PCU_Comm_Begin2(PCUHandle h);
int PCU_Comm_Pack2(PCUHandle h, int to_rank, const void* data, size_t size);
#define PCU_COMM_PACK2(handle, to_rank,object)\
PCU_Comm_Pack2(handle, to_rank,&(object),sizeof(object))
int PCU_Comm_Send2(PCUHandle h);
bool PCU_Comm_Receive2(PCUHandle h);
bool PCU_Comm_Listen2(PCUHandle h);
int PCU_Comm_Sender2(PCUHandle h);
bool PCU_Comm_Unpacked2(PCUHandle h);
int PCU_Comm_Unpack2(PCUHandle h, void* data, size_t size);
#define PCU_COMM_UNPACK2(handle, object)\
PCU_Comm_Unpack2(handle, &(object),sizeof(object))

/*turns deterministic ordering for the
  above API on/off*/
void PCU_Comm_Order2(PCUHandle h, bool on);

/*collective operations*/
void PCU_Barrier2(PCUHandle h);
void PCU_Add_Doubles2(PCUHandle h, double* p, size_t n);
double PCU_Add_Double2(PCUHandle h, double x);
void PCU_Min_Doubles2(PCUHandle h, double* p, size_t n);
double PCU_Min_Double2(PCUHandle h, double x);
void PCU_Max_Doubles2(PCUHandle h, double* p, size_t n);
double PCU_Max_Double2(PCUHandle h, double x);
void PCU_Add_Ints2(PCUHandle h, int* p, size_t n);
int PCU_Add_Int2(PCUHandle h, int x);
void PCU_Add_Longs2(PCUHandle h, long* p, size_t n);
long PCU_Add_Long2(PCUHandle h, long x);
void PCU_Exscan_Ints2(PCUHandle h, int* p, size_t n);
int PCU_Exscan_Int2(PCUHandle h, int x);
void PCU_Exscan_Longs2(PCUHandle h, long* p, size_t n);
long PCU_Exscan_Long2(PCUHandle h, long x);
void PCU_Add_SizeTs2(PCUHandle h, size_t* p, size_t n);
size_t PCU_Add_SizeT2(PCUHandle h, size_t x);
void PCU_Min_SizeTs2(PCUHandle h, size_t* p, size_t n);
size_t PCU_Min_SizeT2(PCUHandle h, size_t x);
void PCU_Max_SizeTs2(PCUHandle h, size_t* p, size_t n);
size_t PCU_Max_SizeT2(PCUHandle h, size_t x);
void PCU_Min_Ints2(PCUHandle h, int* p, size_t n);
int PCU_Min_Int2(PCUHandle h, int x);
void PCU_Max_Ints2(PCUHandle h, int* p, size_t n);
int PCU_Max_Int2(PCUHandle h, int x);
void PCU_Max_Longs2(PCUHandle h, long* p, size_t n);
long PCU_Max_Long2(PCUHandle h, long x);
int PCU_Or2(PCUHandle h, int c);
int PCU_And2(PCUHandle h, int c);

/*process-level self/peers (mpi wrappers)*/
int PCU_Proc_Self2(PCUHandle h);
int PCU_Proc_Peers2(PCUHandle h);

/*IPComMan replacement API*/
int PCU_Comm_Write2(PCUHandle h, int to_rank, const void* data, size_t size);
#define PCU_COMM_WRITE2(handle,to,data) \
PCU_Comm_Write2(handle, to,&(data),sizeof(data))
bool PCU_Comm_Read2(PCUHandle h, int* from_rank, void** data, size_t* size);

/*Debug file I/O API*/
void PCU_Debug_Open2(PCUHandle h);

void PCU_Debug_Print2(PCUHandle h, const char* format, ...) PCU_FORMAT_ATTRIBUTE(2,3);
/*lesser-used APIs*/
bool PCU_Comm_Initialized2(PCUHandle h);
int PCU_Comm_Packed2(PCUHandle h ,int to_rank, size_t* size);
int PCU_Comm_From2(PCUHandle h, int* from_rank);
int PCU_Comm_Received2(PCUHandle h, size_t* size);
void* PCU_Comm_Extract2(PCUHandle h, size_t size);
int PCU_Comm_Rank2(PCUHandle h, int* rank);
int PCU_Comm_Size2(PCUHandle h, int* size);

/*deprecated method enum*/

/*special MPI_Comm replacement API*/
void PCU_Switch_Comm2(PCUHandle h, MPI_Comm new_comm);
MPI_Comm PCU_Get_Comm2(PCUHandle h);

/*stack trace helpers using GNU/Linux*/
void PCU_Protect2(void);

/*MPI_Wtime() equivalent*/
double PCU_Time2(void);

/*Memory usage*/
double PCU_GetMem2(void);

/*Access global variable*/




#ifdef __cplusplus
} /* extern "C" */
#endif

#endif