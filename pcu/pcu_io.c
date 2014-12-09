/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_io.h"
#include "pcu_common.h"
#include "pcu_memory.h"
#include "pcu_mpi.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>

typedef struct pcu_file {
  FILE * f;
  pcu_buffer buf;
  off_t pos;
  bool write;
  bool compress;
} pcu_file;

#ifdef PCU_BZIP
#include <bzlib.h>

static void open_compressed(pcu_file* pf)
{
  off_t file_size;
  fseek(pf->f, 0, SEEK_END);
  file_size = ftello(pf->f);
  rewind(pf->f);

  void* buf = malloc(file_size);
  fread(buf, 1, file_size, pf->f);

  unsigned int len = file_size;
  pcu_make_buffer(&pf->buf);
  pcu_push_buffer(&pf->buf, len);

  int rc;
  while (1)
  {
    rc = BZ2_bzBuffToBuffDecompress (pf->buf.start, &len,
        buf, file_size,
        0, 0);
    if (rc == BZ_OUTBUFF_FULL)
    {
      len += file_size;
      pcu_push_buffer (&pf->buf, file_size);
    }
    else
      break;
  }

  pcu_resize_buffer (&pf->buf, (size_t) len);
  free (buf);

}

static void close_compressed(pcu_file* pf)
{
  void * buf;
  unsigned int len;
  len = pf->buf.size;
  buf = malloc (len);

  /* May change the last value (workfactor) to get different results.*/
  BZ2_bzBuffToBuffCompress (buf, &len,
      pf->buf.start, len,
      1, 0, 0);
  fwrite (buf, 1, len, pf->f);
  free (buf);
}
#else
static void open_compressed(pcu_file* pf)
{
  (void)pf;
  pcu_fail("recompile with bzip2 support");
}

static void close_compressed(pcu_file* pf)
{
  (void)pf;
  pcu_fail("recompile with bzip2 support");
}
#endif

pcu_file* pcu_fopen(const char* name, bool write, bool compress)
{
  pcu_file* pf = (pcu_file*) malloc(sizeof(pcu_file));

  if (write)
    pf->f = fopen(name,"w");
  else
    pf->f = fopen(name,"r");
  if (!pf->f)
  {
    perror("pcu_fopen");
    pcu_fail("fopen failed");
  }

  if (compress)
  {
    if (write)
      pcu_make_buffer(&pf->buf);
    else
      open_compressed(pf);
  }

  pf->compress = compress;
  pf->write = write;
  pf->pos = 0;

  return pf;
}

void pcu_fclose(pcu_file* pf)
{
  // Flush the buffer when writing.
  if (pf->compress && pf->write)
    close_compressed(pf);
  if (pf->compress)
    pcu_free_buffer(&pf->buf);
  fclose(pf->f);
  free(pf);
}

void pcu_fwrite(void const* p, size_t size, size_t nmemb, pcu_file * f)
{
  if (!f->write)
    pcu_fail("file not opened for writing.");

  if (f->compress)
  {
    void* write_me = pcu_push_buffer (&f->buf, size * nmemb);
    memcpy(write_me, p, size * nmemb);
  }
  else
  {
    size_t r = fwrite(p,size,nmemb,f->f);
    if (r != nmemb)
      pcu_fail("fwrite failed");
  }
}

void pcu_fread(void* p, size_t size, size_t nmemb, pcu_file * f)
{
  if (f->write)
    pcu_fail("file not opened for reading.");

  if (f->compress)
  {
    void * read_me = f->buf.start + f->pos;
    memcpy(p, read_me, size * nmemb);
    f->pos += size * nmemb;
  }
  else
  {
    size_t r = fread(p,size,nmemb,f->f);
    if (r != nmemb)
      pcu_fail("fread failed");
  }
}

void pcu_read(pcu_file* f, char* p, size_t n)
{
  pcu_fread(p,1,n,f);
}

void pcu_write(pcu_file* f, const char* p, size_t n)
{
  pcu_fwrite(p,1,n,f);
}

static const uint16_t pcu_endian_value = 1;
#define PCU_ENDIANNESS ((*((uint8_t*)(&pcu_endian_value)))==1)
#define PCU_BIG_ENDIAN 0
#define PCU_ENCODED_ENDIAN PCU_BIG_ENDIAN //consistent with network byte order

static void pcu_swap_16(uint16_t* p)
{
  uint8_t* p2 = (uint8_t*)p;
  uint8_t temp = p2[0];
  p2[0] = p2[1];
  p2[1] = temp;
}

static void pcu_swap_32(uint32_t* p)
{
  uint16_t* p2 = (uint16_t*)p;
  uint16_t temp = p2[0];
  p2[0] = p2[1];
  p2[1] = temp;
  pcu_swap_16(p2);
  pcu_swap_16(p2+1);
}

static void pcu_swap_64(uint32_t* p)
{
  uint32_t temp = p[0];
  p[0] = p[1];
  p[1] = temp;
  pcu_swap_32(p);
  pcu_swap_32(p+1);
}

void pcu_swap_unsigneds(unsigned* p, size_t n)
{
  assert(sizeof(unsigned)==4);
  for (size_t i=0; i < n; ++i)
    pcu_swap_32(p++);
}

void pcu_swap_doubles(double* p, size_t n)
{
  assert(sizeof(double)==8);
  for (size_t i=0; i < n; ++i)
    pcu_swap_64((uint32_t*)(p++));
}

void pcu_write_unsigneds(pcu_file* f, unsigned* p, size_t n)
{
  unsigned* tmp;
  if (PCU_ENDIANNESS != PCU_ENCODED_ENDIAN) {
    tmp = malloc(n * sizeof(unsigned));
    memcpy(tmp, p, n * sizeof(unsigned));
    pcu_swap_unsigneds(tmp, n);
    pcu_fwrite(tmp,sizeof(unsigned),n,f);
    free(tmp);
  } else {
    pcu_fwrite(p,sizeof(unsigned),n,f);
  }
}

void pcu_write_doubles(pcu_file* f, double* p, size_t n)
{
  double* tmp;
  if (PCU_ENDIANNESS != PCU_ENCODED_ENDIAN) {
    tmp = malloc(n * sizeof(double));
    memcpy(tmp, p, n * sizeof(double));
    pcu_swap_doubles(tmp, n);
    pcu_fwrite(tmp,sizeof(double),n,f);
    free(tmp);
  } else {
    pcu_fwrite(p,sizeof(double),n,f);
  }
}

void pcu_read_unsigneds(pcu_file* f, unsigned* p, size_t n)
{
  pcu_fread(p,sizeof(unsigned),n,f);
  if (PCU_ENDIANNESS != PCU_ENCODED_ENDIAN)
    pcu_swap_unsigneds(p,n);
}

void pcu_read_doubles(pcu_file* f, double* p, size_t n)
{
  pcu_fread(p,sizeof(double),n,f);
  if (PCU_ENDIANNESS != PCU_ENCODED_ENDIAN)
    pcu_swap_doubles(p,n);
}

void pcu_read_string (pcu_file* f, char ** p)
{
  pcu_buffer buf;
  pcu_make_buffer(&buf);
  char* c;
  do {
    pcu_push_buffer(&buf,1);
    c = buf.start + buf.size - 1;
    pcu_read(f,c,1);
  } while (*c != '\0');
  pcu_resize_buffer(&buf,buf.size);
  *p = buf.start;
}

void pcu_write_string (pcu_file * f, const char * p)
{
  size_t len = strlen (p);
  pcu_write (f, p, len + 1);
}

FILE* pcu_open_parallel(const char* prefix, const char* ext)
{
  //max_rank_chars = strlen("4294967296"), 4294967296 = 2^32 ~= INT_MAX
  static const size_t max_rank_chars = 10;
  size_t path_size = strlen(prefix) + max_rank_chars + strlen(ext) + 1;
  char* path = pcu_malloc(path_size);
  int rank = pcu_mpi_rank();
  snprintf(path,path_size,"%s%d.%s",prefix,rank,ext);
  FILE* file = fopen(path, "w");
  pcu_free(path);
  return file;
}
