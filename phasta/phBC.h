
#ifndef PH_BC_H
#define PH_BC_H

#include <map>
#include <set>
#include <string>

namespace ph {

struct BC
{
  BC();
  ~BC();
  int tag;
  int dim;
  double* values;
  bool operator<(const BC& other) const;
};

struct FieldBCs
{
  int size;
  typedef std::set<BC> Set;
  Set bcs;
};

struct BCs
{
  typedef std::map<std::string, FieldBCs> Map;
  Map fields;
};

void readBCs(const char* filename, BCs& bcs);

}

#endif
