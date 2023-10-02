/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_io.h"
#include "noto_malloc.h"
#include "reel.h"
#include "pcu_mpi.h"
#include "PCU.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "pcu_util.h"
#include <sys/types.h>
#include <limits.h>

#ifdef PCU_BZIP
#include <bzlib.h>
#endif

typedef struct pcu_file {
  FILE* f;
#ifdef PCU_BZIP
  BZFILE* bzf;
#endif
  bool write;
  bool compress;
} pcu_file;

#ifdef PCU_BZIP

static void open_compressed_read(pcu_file* pf)
{
  int bzerror;
  int verbosity = 0;
  int small = 0;
  void* unused = NULL;
  int nUnused = 0;
  pf->bzf = BZ2_bzReadOpen(&bzerror, pf->f, verbosity, small, unused, nUnused);
  if (bzerror != BZ_OK)
    reel_fail("BZ2_bzReadOpen failed with code %d", bzerror);
}

static void open_compressed_write(pcu_file* pf)
{
  int bzerror;
  int verbosity = 0;
  int blockSize100k = 9;
  int workFactor = 30;
  pf->bzf = BZ2_bzWriteOpen(&bzerror, pf->f, blockSize100k, verbosity, workFactor);
  if (bzerror != BZ_OK)
    reel_fail("BZ2_bzWriteOpen failed with code %d", bzerror);
}

static void open_compressed(pcu_file* pf)
{
  if (pf->write)
    open_compressed_write(pf);
  else
    open_compressed_read(pf);
}

static void compressed_read(pcu_file* pf, void* data, size_t size)
{
  int bzerror;
  int len;
  int rv;
  PCU_ALWAYS_ASSERT(size < INT_MAX);
  len = size;
  rv = BZ2_bzRead(&bzerror, pf->bzf, data, len);
  if (bzerror != BZ_OK && bzerror != BZ_STREAM_END)
    reel_fail("BZ2_bzRead failed with code %d", bzerror);
  PCU_ALWAYS_ASSERT(rv == len);
}

static void compressed_write(pcu_file* pf, void const* data, size_t size)
{
  int bzerror;
  int len;
  void* bzip2_is_not_const_correct;
  PCU_ALWAYS_ASSERT(size < INT_MAX);
  len = size;
  bzip2_is_not_const_correct = (void*)data;
  BZ2_bzWrite(&bzerror, pf->bzf, bzip2_is_not_const_correct, len);
  if (bzerror != BZ_OK)
    reel_fail("BZ2_bzWrite failed with code %d", bzerror);
}

static void close_compressed_read(pcu_file* pf)
{
  int bzerror;
  BZ2_bzReadClose(&bzerror, pf->bzf);
  if (bzerror != BZ_OK)
    reel_fail("BZ2_readClose failed with code %d", bzerror);
}

static void close_compressed_write(pcu_file* pf)
{
  int bzerror;
  int abandon = 0;
  unsigned* nbytes_in = NULL;
  unsigned* nbytes_out = NULL;
  BZ2_bzWriteClose(&bzerror, pf->bzf, abandon, nbytes_in, nbytes_out);
  if (bzerror != BZ_OK)
    reel_fail("BZ2_writeClose failed with code %d", bzerror);
}

static void close_compressed(pcu_file* pf)
{
  if (pf->write)
    close_compressed_write(pf);
  else
    close_compressed_read(pf);
}

#else

static void open_compressed(pcu_file* pf)
{
  (void)pf;
  reel_fail("recompile PCU with -DPCU_COMPRESS=ON");
}

static void compressed_read(pcu_file* pf, void* data, size_t size)
{
  (void)pf;
  (void)data;
  (void)size;
  reel_fail("recompile PCU with -DPCU_COMPRESS=ON");
}

static void compressed_write(pcu_file* pf, void const* data, size_t size)
{
  (void)pf;
  (void)data;
  (void)size;
  reel_fail("recompile PCU with -DPCU_COMPRESS=ON");
}

static void close_compressed(pcu_file* pf)
{
  (void)pf;
  reel_fail("recompile PCU with -DPCU_COMPRESS=ON");
}

#endif

/**
 * brief limit the number of ranks that can call fopen simultaneously
 * remark Argonne's GPFS filesystem is failing to open some files when
 *        tens of thousands of ranks simultaneously make the request.
 *        Ideally, the filesystem handles the load and this can be
 *        removed.
 */
FILE* pcu_group_open(const char* path, bool write) {
  FILE* fp = NULL;
  const int rank = PCU_Comm_Self();
  const char* mode = write ? "w" : "r";
  const int group_size = 4096;
  const int q = PCU_Comm_Peers()/group_size;
  const int r = PCU_Comm_Peers()%group_size;
  const int groups = q + ( r > 0 );
  if(!rank && groups > 1) {
    fprintf(stderr,
        "pcu peers %d max group size %d posix groups %d\n",
        PCU_Comm_Peers(), group_size, groups);
  }
  for(int i=0; i<groups; i++) {
    if(rank%groups == i) {
      fp = fopen(path, mode);
      if (!fp)
        reel_fail("Could not find or open file \"%s\"\n", path);
    }
    PCU_Barrier();
  }
  return fp;
}

pcu_file* pcu_fopen(const char* name, bool write, bool compress)
{
  pcu_file* pf = (pcu_file*) malloc(sizeof(pcu_file));
  pf->compress = compress;
  pf->write = write;
  pf->f = pcu_group_open(name, write);
  if (!pf->f) {
    perror("pcu_fopen");
    reel_fail("pcu_fopen couldn't open \"%s\"", name);
  }
  if(compress)
    open_compressed(pf);
  return pf;
}

void pcu_fclose(pcu_file* pf)
{
  if (pf->compress)
    close_compressed(pf);
  fclose(pf->f);
  free(pf);
}

void pcu_fwrite(void const* p, size_t size, size_t nmemb, pcu_file * f)
{
  if (!f->write)
    reel_fail("pcu_fwrite: file not opened for writing.");
  if (f->compress) {
    compressed_write(f, p, size * nmemb);
  } else {
    if (nmemb != fwrite(p, size, nmemb, f->f))
      reel_fail("fwrite(%p, %lu, %lu, %p) failed", p, size, nmemb, (void*) f->f);
  }
}

void pcu_fread(void* p, size_t size, size_t nmemb, pcu_file * f)
{
  if (f->write)
    reel_fail("pcu_fread: file not opened for reading.");
  if (f->compress) {
    compressed_read(f, p, size * nmemb);
  } else {
    if (nmemb != fread(p, size, nmemb, f->f))
      reel_fail("fread(%p, %lu, %lu, %p) failed", p, size, nmemb, (void*) f->f);
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
  PCU_ALWAYS_ASSERT(sizeof(unsigned)==4);
  for (size_t i=0; i < n; ++i)
    pcu_swap_32(p++);
}

void pcu_swap_doubles(double* p, size_t n)
{
  PCU_ALWAYS_ASSERT(sizeof(double)==8);
  for (size_t i=0; i < n; ++i)
    pcu_swap_64((uint32_t*)(p++));
}

void pcu_write_unsigneds(pcu_file* f, unsigned* p, size_t n)
{
  unsigned* tmp;
  if (n)
    PCU_ALWAYS_ASSERT(p != 0);
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
  if (n)
    PCU_ALWAYS_ASSERT(p != 0);
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
  char* path = noto_malloc(path_size);
  int rank = pcu_mpi_rank();
  snprintf(path,path_size,"%s%d.%s",prefix,rank,ext);
  FILE* file = fopen(path, "w");
  noto_free(path);
  return file;
}
