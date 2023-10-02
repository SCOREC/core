#include <lionPrint.h>
#include <cstdio>
#include <cassert>
#include <cstdlib>

int main(int argc, char** argv) {
  if( argc != 2 ) {
    fprintf(stderr, "Usage: %s <verbosity level>\n", argv[0]);
    return 0;
  }
  int lvl=atoi(argv[1]);
  lion_set_verbosity(lvl);
  lion_oprint(1,"message lvl 1\n");
  lion_oprint(2,"message lvl 2\n");
  lion_eprint(3,"message lvl 3\n");
  lion_set_verbosity(0);
  int ret = lion_oprint(1,"if you see this message"
                  "there is a problem and the test"
                  "has failed\n");
  if( ret ) return 1;
  FILE* f = fopen("outputGoesHere.txt","w");
  lion_set_stdout(f);
  lion_set_stderr(f);
  lion_set_verbosity(lvl);
  lion_oprint(1,"out message lvl 1\n");
  lion_eprint(1,"err message lvl 1\n");
  fclose(f);
  return 0;
}
