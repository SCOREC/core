#ifndef PARMA_ASSOCIATIVE_H
#define PARMA_ASSOCIATIVE_H
#include <map>
#include <sstream>
#include <pcu_util.h>

namespace parma {
  template <class T> class Associative {
    typedef std::map<int, T> Container;
    public:
      Associative() {
        iteratorActive = false;
      }
      void begin() {
        PCU_ALWAYS_ASSERT(!iteratorActive);
        iteratorActive = true;
        cItr = c.begin();
      }
      typedef std::pair<const int, T> Item;
      const Item* iterate() {
        PCU_ALWAYS_ASSERT(iteratorActive);
        if( cItr == c.end() ) 
          return NULL;
        else
          return &(*cItr++); // there is no spoon ... and this is likely crap
      }
      void end() {
        PCU_ALWAYS_ASSERT(iteratorActive);
        iteratorActive = false;
      }
      T get(int key) {
        return c[key];
      }
      void set(int key, T value) {
        c[key] = value;
      }
      bool has(int key) {
        return (c.count(key) != 0);
      }
      size_t size() {
        return c.size();
      }
      std::string print(const char* key) {
        std::stringstream s;
        s << key << " ";
        const Item* i;
        begin();
        while( (i = iterate()) ) 
          s << i->first << " " << i->second << " ";
        end();
        return s.str();
      }
    protected:
      Container c;
    private:
      typename Container::iterator cItr;
      bool iteratorActive;
  };
}
#endif
