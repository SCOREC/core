#ifndef PH_INPUT_H
#define PH_INPUT_H

#include <string>

struct RStream;

namespace ph {

class Input
{
  public:
    Input();
    void load(const char* filename);
    int globalP;
    int timeStepNumber;
    int ensa_dof;
    std::string restartFileName;
    std::string attributeFileName;
    std::string meshFileName;
    std::string outMeshFileName;
    std::string modelFileName;
    std::string outputFormat;
    std::string partitionMethod;
    std::string preAdaptBalanceMethod;
    std::string prePhastaBalanceMethod;
    int adaptFlag;
    int rRead;
    int rStart;
    int adaptStrategy;
    double adaptErrorThreshold;
    std::string adaptErrorFieldName;
    int adaptErrorFieldIndex;
    int periodic;
    int prCD;
    int timing;
    int internalBCNodes;
    int writeDebugFiles;
    int phastaIO;
    int splitFactor;
    int solutionMigration;
    int displacementMigration;
    int isReorder;
    int tetrahedronize;
    int localPtn;
    int recursiveUR;
    int dwalMigration;
    int buildMapping;
    int elementsPerMigration;
    int threaded;
    int initBubbles;
    int formElementGraph;
    int snap;
    int splitAllLayerEdges;
    int filterMatches;
    int axisymmetry;
    int formEdges;
    double elementImbalance;
    double vertexImbalance;
    FILE* (*openfile_read)(Input& in, const char* path);
    RStream* rs;
};

int countNaturalBCs(Input& in);
int countEssentialBCs(Input& in);
int countScalarBCs(Input& in);

}

#endif


