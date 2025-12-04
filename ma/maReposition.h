#ifndef MA_REPOSITION_H
#define MA_REPOSITION_H

#include "maMesh.h"
#include "ma.h"

namespace ma {

class RepositionVertex
{
    public:
    RepositionVertex(Adapt* a);
    bool move(Entity* vertex, Vector target);
    void cancel(Entity* vertex);
    apf::Up& getInvalid();
    double getWorstQuality();

    private:
    Adapt* adapt;
    Mesh* mesh;
    Entity* vertex;
    Vector prevPosition;
    Upward adjacentElements;
    apf::Up invalid;
    double worstQuality;

    void findInvalid();
};

bool repositionVertex(Mesh* m, Entity* v,
    int max_iters, double initial_speed);

}

#endif
