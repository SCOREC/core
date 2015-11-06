#ifndef SAMSZ_H
#define SAMSZ_h
namespace apf { 
  class Field; 
  class Mesh;
}
namespace samSz {
  apf::Field* isoSize(apf::Mesh* m);
}
#endif
