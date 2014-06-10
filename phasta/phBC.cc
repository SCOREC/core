#include "phBC.h"
#include <apf.h>
#include <fstream>
#include <sstream>

namespace ph {

BC::BC()
{
  values = 0;
}

BC::~BC()
{
  delete [] values;
}

bool BC::operator<(const BC& other) const
{
  if (dim != other.dim)
    return dim < other.dim;
  return tag < other.tag;
}

static struct { const char* name; int size; } builtin_sizes[4] =
{{"initial velocity", 3}
,{"comp1", 4}
,{"comp3", 4}
,{"traction vector", 3}
};

static int getSize(std::string const& name)
{
  for (int i = 0; i < 4; ++i)
    if (name == builtin_sizes[i].name)
      return builtin_sizes[i].size;
  return 1;
}

static void readBC(std::string const& line, BCs& bcs)
{
  std::stringstream ss(line);
  std::string name;
  std::getline(ss, name, ':');
  if (!bcs.fields.count(name)) {
    FieldBCs fbcs;
    fbcs.size = getSize(name);
    bcs.fields[name] = fbcs;
  }
  FieldBCs& fbcs = bcs.fields[name];
  BC bc;
  ss >> bc.tag >> bc.dim;
  bc.values = new double[fbcs.size];
  for (int i = 0; i < fbcs.size; ++i)
    ss >> bc.values[i];
  fbcs.bcs.insert(bc);
  bc.values = 0; //ownership of pointer transferred, prevent delete from here
}

void readBCs(const char* filename, BCs& bcs)
{
  std::ifstream file(filename);
  std::string line;
  while (std::getline(file, line, '\n')) {
    if (line[0] == '#')
      continue;
    readBC(line, bcs);
  }
}

}
