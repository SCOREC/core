/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_IO_H
#define PCU_IO_H


#ifdef __cplusplus
#include <cstdio>
extern "C" {
#else
#include <stdio.h>
#include <stdbool.h>
#endif

struct pcu_file;
typedef struct PCU_t PCU_t;

struct pcu_file* pcu_fopen(PCU_t h, const char* path, bool write, bool compress);
void pcu_fclose (struct pcu_file * pf);
void pcu_read(struct pcu_file* f, char* p, size_t n);
void pcu_write(struct pcu_file* f, const char* p, size_t n);
void pcu_read_unsigneds(struct pcu_file* f, unsigned* p, size_t n);
#define PCU_READ_UNSIGNED(f,p) pcu_read_unsigneds(f,&(p),1);
void pcu_write_unsigneds(struct pcu_file* f, unsigned* p, size_t n);
#define PCU_WRITE_UNSIGNED(f,p) pcu_write_unsigneds(f,&(p),1);
void pcu_read_doubles(struct pcu_file* f, double* p, size_t n);
void pcu_write_doubles(struct pcu_file* f, double* p, size_t n);
void pcu_read_string(struct pcu_file* f, char** p);
void pcu_write_string(struct pcu_file* f, const char* p);

FILE* pcu_open_parallel(PCU_t h, const char* prefix, const char* ext);
FILE* pcu_group_open(PCU_t h, const char* path, bool write);

void pcu_swap_doubles(double* p, size_t n);
void pcu_swap_unsigneds(unsigned* p, size_t n);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
