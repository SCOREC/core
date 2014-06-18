#ifndef PARMA_BASE_H
#define PARMA_BASE_H
#include <map>

namespace parma {
  template <class T> class Associative {
    typedef std::map<int, T> Container;
    public:
      Associative();
      void begin();
      typedef std::pair<const int, T> Item;
      const Item* iterate();
      void end();
      T get(int key);
      void set(int key, T value);
      bool has(int key);
      void print(const char* key);
    protected:
      Container c;
    private:
      typename Container::iterator cItr;
      bool iteratorActive;
  };
  

};
#endif
