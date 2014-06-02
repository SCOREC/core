#include <phBC.h>
#include <apf.h>

int main(int argc, char** argv)
{
  ph::BCs bcs;
  ph::readBCs(argv[1], bcs);
  APF_ITERATE(ph::BCs::Map, bcs.fields, it) {
    ph::FieldBCs& fbc = it->second;
    APF_ITERATE(ph::FieldBCs::Set, fbc.bcs, it2) {
      std::cout << it->first << ": ";
      std::cout << it2->tag << ' ';
      std::cout << it2->dim << ' ';
      for (int i = 0; i < fbc.size; ++i)
        std::cout << it2->values[i] << ' ';
      std::cout << '\n';
    }
  }
}
