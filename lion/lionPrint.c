#include "lionPrint.h"
#include <assert.h>
#include <stdarg.h>

int lion_verbosity_level = 0;
FILE* lion_stdout = NULL;
FILE* lion_stderr = NULL;

int lion_get_verbosity() {
  return lion_verbosity_level;
}

void lion_set_verbosity(int lvl) {
  assert(lvl >= 0);
  lion_verbosity_level = lvl;
}

void lion_set_stdout(FILE* out) {
  assert(out != NULL);
  lion_stdout = out;
}

void lion_set_stderr(FILE* err) {
  assert(err != NULL);
  lion_stderr = err;
}

#define print(lvl,dest,fmt) \
  va_list ap; \
  va_start(ap,fmt); \
  int written = 0; \
  if(lion_verbosity_level >= lvl) \
    written = vfprintf(dest,fmt,ap); \
  va_end(ap); \
  return written;

#define setstream(stream,defaultStream) \
  static int calls = 0; \
  if( !calls && !stream) \
    stream = defaultStream; \
  calls++; \

int lion_oprint(int lvl, char const* fmt, ...) {
  setstream(lion_stdout,stdout);
  assert(lvl >= 0);
  print(lvl,lion_stdout,fmt);
}

int lion_eprint(int lvl, char const* fmt, ...) {
  setstream(lion_stderr,stderr);
  assert(lvl >= 0);
  print(lvl,lion_stderr,fmt);
}

#define vprint(lvl,dest,fmt,ap) \
  int written = 0; \
  if(lion_verbosity_level >= lvl) \
    written = vfprintf(dest,fmt,ap); \
  return written;

int lion_voprint(int lvl, char const* fmt, va_list ap) {
  setstream(lion_stdout,stdout);
  assert(lvl >= 0);
  vprint(lvl,lion_stdout,fmt,ap);
}

int lion_veprint(int lvl, char const* fmt, va_list ap) {
  setstream(lion_stderr,stderr);
  assert(lvl >= 0);
  vprint(lvl,lion_stderr,fmt,ap);
}
