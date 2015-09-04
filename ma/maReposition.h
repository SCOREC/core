#ifndef MA_REPOSITION_H
#define MA_REPOSITION_H

#include "maMesh.h"

namespace ma {

bool repositionVertex(Mesh* m, Entity* v,
    int max_iters, double initial_speed);

}

#endif
