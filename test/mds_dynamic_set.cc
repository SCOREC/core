#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif
#include <mds.h>
#ifdef __cplusplus
}
#endif

int main() {
  struct mds_set s;
  mds_init_set(&s);
  mds_expand_set(&s, 1024);
  if (s.cap < 1024) {
    std::cerr << "ERROR: s->cap < 1024." << std::endl;
  }
  mds_expand_set(&s, 2934);
  if (s.cap < 2934) {
    std::cerr << "ERROR: s->cap < 1024." << std::endl;
  }
  mds_destroy_set(&s);
  return 0;
}
