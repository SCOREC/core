#ifndef SAMSZ_H
#define SAMSZ_H
namespace apf { 
  class Field; 
  class Mesh;
}
namespace samSz {
  apf::Field* isoSize(apf::Mesh* m);
}
#endif
