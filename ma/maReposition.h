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
    bool moveToImproveQuality(Entity* vertex);
    void moveToImproveShortEdges(Entity* vertex);
    bool improveQuality(Entity* tet);
    void cancel(Entity* vertex);
    apf::Up& getInvalid();

    private:
    Adapt* adapt;
    Mesh* mesh;
    Entity* vertex;
    Vector prevPosition;
    Upward adjacentElements;
    apf::Up invalid;
    apf::Up adjEdges;
    double worstQuality;
    double startingQuality;

    void findInvalid();
    void clearAdjCache();
    Vector modelCenter();
    void init(Entity* vertex);
    double findWorstShape(Vector position);
    double findShortestEdge(Vector position);
};

bool repositionVertex(Mesh* m, Entity* v,
    int max_iters, double initial_speed);

}

#endif
