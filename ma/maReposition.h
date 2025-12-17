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
    void moveToHighestQuality(Entity* vertex);
    void cancel(Entity* vertex);
    apf::Up& getInvalid();

    private:
    Adapt* adapt;
    Mesh* mesh;
    Entity* vertex;
    Vector prevPosition;
    Upward adjacentElements;
    apf::Up invalid;
    double worstQuality;
    std::vector<double> oldCache;

    void findInvalid();
    void storeOldCache();
    Vector cavityCenter();
    double findWorstShape(Vector position);
};

bool repositionVertex(Mesh* m, Entity* v,
    int max_iters, double initial_speed);

}

#endif
