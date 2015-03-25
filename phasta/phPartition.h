#ifndef PH_PARTITION
#define PH_PARTITION

namespace apf {
class Mesh2;
class Migration;
}

namespace ph {

class Input;

void split(Input& in, apf::Mesh2* m, void (*runAfter)(apf::Mesh2*));
apf::Migration* split(Input& in, apf::Mesh2* m);
void balance(apf::Mesh2* m);

}

#endif
