/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PUMI_COMM_H
#define PUMI_COMM_H

#define PUMI_COMM_FAIL -42

#ifdef __cplusplus
#include <cstddef>
extern "C" {
#else
#include <stddef.h>
#include <stdbool.h>
#endif

#ifdef __cplusplus
enum { PUMI_COMM_GLOBAL, PUMI_COMM_LOCAL };
#else
enum { PUMI_COMM_GLOBAL, PUMI_COMM_LOCAL };
#endif

typedef void* (*PUMI_Thrd_Func)(void*);

int PUMI_Thrd_Run(int nthreads, PUMI_Thrd_Func function, void** in_out);
int PUMI_Thrd_Rank(int* rank);
int PUMI_Thrd_Size(int* size);

int PUMI_Comm_Init(void);
int PUMI_Comm_Rank(int* rank);
int PUMI_Comm_Size(int* size);
int PCU_Comm_Self(void);
#define PUMI_Comm_Self PCU_Comm_Self
int PCU_Comm_Peers(void);
#define PUMI_Comm_Peers PCU_Comm_Peers
int PUMI_Comm_Start(int method);
void PCU_Comm_Begin(void);
#define PUMI_Comm_Begin PCU_Comm_Begin
int PUMI_Comm_Pack(int to_rank, const void* data, size_t size);
#define PUMI_COMM_PACK(to_rank,object)\
PUMI_Comm_Pack(to_rank,&(object),sizeof(object))
int PUMI_Comm_Packed(int to_rank, size_t* size);
int PUMI_Comm_Write(int to_rank, const void* data, size_t size);
#define PUMI_COMM_WRITE(to,data) \
PUMI_Comm_Write(to,&(data),sizeof(data))
int PUMI_Comm_Send(void);
int PUMI_Comm_Receive(bool* done);
bool PCU_Comm_Listen(void);
#define PUMI_Comm_Listen PCU_Comm_Listen
bool PUMI_Comm_Read(int* from_rank, void** data, size_t* size);
int PUMI_Comm_From(int* from_rank);
int PCU_Comm_Sender(void);
#define PUMI_Comm_Sender PCU_Comm_Sender
int PUMI_Comm_Received(size_t* size);
int PUMI_Comm_Unpack(void* data, size_t size);
#define PUMI_COMM_UNPACK(object)\
PUMI_Comm_Unpack(&(object),sizeof(object))
bool PUMI_Comm_Unpacked(void);
void* PCU_Comm_Extract(size_t size);
#define PUMI_Comm_Extract PCU_Comm_Extract
int PUMI_Comm_Free(void);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
