#ifndef PCU_C_H
#define PCU_C_H
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

typedef struct PCU_t PCU_t;

int PCU_Comm_Init(PCU_t* h);
int PCU_Comm_Free(PCU_t* h);

/*rank/size functions*/
int PCU_Comm_Self(PCU_t h);
int PCU_Comm_Peers(PCU_t h);

/*recommended message passing API*/
void PCU_Comm_Begin(PCU_t h);
int PCU_Comm_Pack(PCU_t h, int to_rank, const void* data, size_t size);
#define PCU_COMM_PACK(handle, to_rank,object)\
PCU_Comm_Pack(handle, to_rank,&(object),sizeof(object))
int PCU_Comm_Send(PCU_t h);
bool PCU_Comm_Receive(PCU_t h);
bool PCU_Comm_Listen(PCU_t h);
int PCU_Comm_Sender(PCU_t h);
bool PCU_Comm_Unpacked(PCU_t h);
int PCU_Comm_Unpack(PCU_t h, void* data, size_t size);
#define PCU_COMM_UNPACK(handle, object)\
PCU_Comm_Unpack(handle, &(object),sizeof(object))

/*turns deterministic ordering for the
  above API on/off*/
void PCU_Comm_Order(PCU_t h, bool on);

/*collective operations*/
void PCU_Barrier(PCU_t h);
void PCU_Add_Doubles(PCU_t h, double* p, size_t n);
double PCU_Add_Double(PCU_t h, double x);
void PCU_Min_Doubles(PCU_t h, double* p, size_t n);
double PCU_Min_Double(PCU_t h, double x);
void PCU_Max_Doubles(PCU_t h, double* p, size_t n);
double PCU_Max_Double(PCU_t h, double x);
void PCU_Add_Ints(PCU_t h, int* p, size_t n);
int PCU_Add_Int(PCU_t h, int x);
void PCU_Add_Longs(PCU_t h, long* p, size_t n);
long PCU_Add_Long(PCU_t h, long x);
void PCU_Exscan_Ints(PCU_t h, int* p, size_t n);
int PCU_Exscan_Int(PCU_t h, int x);
void PCU_Exscan_Longs(PCU_t h, long* p, size_t n);
long PCU_Exscan_Long(PCU_t h, long x);
void PCU_Add_SizeTs(PCU_t h, size_t* p, size_t n);
size_t PCU_Add_SizeT(PCU_t h, size_t x);
void PCU_Min_SizeTs(PCU_t h, size_t* p, size_t n);
size_t PCU_Min_SizeT(PCU_t h, size_t x);
void PCU_Max_SizeTs(PCU_t h, size_t* p, size_t n);
size_t PCU_Max_SizeT(PCU_t h, size_t x);
void PCU_Min_Ints(PCU_t h, int* p, size_t n);
int PCU_Min_Int(PCU_t h, int x);
void PCU_Max_Ints(PCU_t h, int* p, size_t n);
int PCU_Max_Int(PCU_t h, int x);
void PCU_Max_Longs(PCU_t h, long* p, size_t n);
long PCU_Max_Long(PCU_t h, long x);
int PCU_Or(PCU_t h, int c);
int PCU_And(PCU_t h, int c);

/*process-level self/peers (mpi wrappers)*/
int PCU_Proc_Self(PCU_t h);
int PCU_Proc_Peers(PCU_t h);

/*IPComMan replacement API*/
int PCU_Comm_Write(PCU_t h, int to_rank, const void* data, size_t size);
#define PCU_COMM_WRITE(handle,to,data) \
PCU_Comm_Write(handle, to,&(data),sizeof(data))
bool PCU_Comm_Read(PCU_t h, int* from_rank, void** data, size_t* size);

/*Debug file I/O API*/
void PCU_Debug_Open(PCU_t h);

void PCU_Debug_Print(PCU_t h, const char* format, ...) PCU_FORMAT_ATTRIBUTE(2,3);
/*lesser-used APIs*/
bool PCU_Comm_Initialized(PCU_t h);
int PCU_Comm_Packed(PCU_t h ,int to_rank, size_t* size);
int PCU_Comm_From(PCU_t h, int* from_rank);
int PCU_Comm_Received(PCU_t h, size_t* size);
void* PCU_Comm_Extract(PCU_t h, size_t size);
int PCU_Comm_Rank(PCU_t h, int* rank);
int PCU_Comm_Size(PCU_t h, int* size);

/*deprecated method enum*/

/*special MPI_Comm replacement API*/
void PCU_Switch_Comm(PCU_t h, MPI_Comm new_comm);
MPI_Comm PCU_Get_Comm(PCU_t h);

/*stack trace helpers using GNU/Linux*/
void PCU_Protect(void);

/*MPI_Wtime() equivalent*/
double PCU_Time(void);

/*Memory usage*/
double PCU_GetMem(void);

/*Access global variable*/
PCU_t PCU_Get_Global_Handle(void);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif