#ifndef PH_PARTITION
#define PH_PARTITION

namespace apf {
class Mesh2;
class Migration;
}

namespace ph {

class Input;

apf::Migration* split(Input& in, apf::Mesh2* m);
void balance(Input& in, apf::Mesh2* m);

}

#endif
