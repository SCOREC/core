#ifndef PH_SPLIT
#define PH_SPLIT

namespace apf {
class Mesh2;
}

namespace ph {

class Input;

void split(Input& in, apf::Mesh2* m, void (*runAfter)(apf::Mesh2*));

}

#endif
