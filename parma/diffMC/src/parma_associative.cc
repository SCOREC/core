#include "parma_associative.h"
#include <assert.h>

namespace parma {
  template <class T> Associative<T>::Associative() {
    iteratorActive = false;
  }
  void Associative::begin() {
    assert(!iteratorActive);
    iteratorActive = true;
    cItr = c.begin();
  }
  const Associative::Item* Associative::iterate() {
    assert(iteratorActive);
    if( cItr == c.end() ) 
      return NULL;
    else
      return &(*cItr++); // there is no spoon ... and this is likely crap
  }
  void Associative::end() {
    assert(iteratorActive);
    iteratorActive = false;
  }
  T Associative::get(int key) {
    return c[key];
  }
  void Associative::set(int key, T value) {
    c[key] = value;
  }
  bool Associative::has(int key) {
    return (c.count(key) != 0);
  }
  void Associative::print(const char* key) {
    std::stringstream s;
    s << key << " ";
    const Item* i;
    begin();
    while( (i = iterate()) ) 
      s << i->first << " " << i->second << " ";
    end();
    std::string str = s.str();
    fprintf(stdout, "%s\n", str.c_str());
  }
}
