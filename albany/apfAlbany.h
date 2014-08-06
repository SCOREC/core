#ifndef APF_ALBANY_H

namespace apf {
class Mesh2;
}

namespace alb {

void shrinkPartition(apf::Mesh2* m, int factor,
    void (*runAfter)(apf::Mesh2* m));

}

#endif
