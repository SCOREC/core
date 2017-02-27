#ifndef PH_INPUT_H
#define PH_INPUT_H

/** \file phInput.
 * h
    \brief The Chef interface for execution control
    \details The variables defined here should be placed in a file named
            'adapt.inp'.  Each variable should be on its own line and
            followed by a space and then the value assigned to the
            variable. Blank lines and lines starting with '#' are ignored.
            To add a new variable edit this file and phInput.cc .
*/

#include <string>

struct RStream;

namespace ph {

/** \brief User configuration for Chef execution */
class Input
{
  public:
    Input();
    void load(const char* filename);
    int timeStepNumber;
    /** \brief this corresponds to the number of degrees of
      freedom in the solution field of the output restart file.
      Note that it should correspond to the number of initial
      conditions specified in the spj file if the solution is
      built from scratch. When the solution is migrated from
      existing restart files, it should also correspond to the
      number of dof in the existing solution field. Set it to 5
      for single phase flow with no turbulence model and 7 for
      two phase flows with level set scalars. */
    int ensa_dof;
    int ensa_melas_dof;
    /** \brief  path to the restart files
        \details this will be read in when solution migration is activated.
      The path should be written as 'N-procs_case/restart' where N is the
      number of processes. The phasta reader will then add the time step
      stamp to the name of this restartFileName variable, as well as the file
      number. */
    std::string restartFileName;
    /** \brief path to the spj or smd file containing the boundary
       and initial conditions*/
    std::string attributeFileName;
    /** \brief path to the directory that includes the input mesh
        \details the path to the SCOREC MDS mesh must end with a '/' if it is a
       directory containing multiple '<partid>.smb' files. This path
       can also be prepended by "bz2:" to tell the mesh file reader that the
       files have been compressed. */
    std::string meshFileName;
    /** \brief output mesh file name, see meshFileName */
    std::string outMeshFileName;
    /** \brief path to the geometric model
        \details the SCOREC discrete (.dmg) is supported if core was built
       without Simmetrix enabled. If Simmetrix is enabled Parasolid (.x_t),
       ACIS (.sat), and Simmetrix GeomSim (.smd) models are supported */
    std::string modelFileName;
    std::string outputFormat;
    std::string partitionMethod;
    std::string preAdaptBalanceMethod;
    std::string midAdaptBalanceMethod;
    std::string postAdaptBalanceMethod;
    std::string prePhastaBalanceMethod;
    int adaptFlag;
    int rRead;
    int rStart;
    int adaptStrategy;
    double validQuality;
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
    int useAttachedFields;
    int displacementMigration;
    int isReorder;
/** \brief tetrahedronize a mixed mesh if set to 1. */
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
    int transferParametric;
    int splitAllLayerEdges;
/** \brief filter out a subset of 3-way periodic matches.
   \ it also fileter out DG ineterface matches. */
    int filterMatches;
    int axisymmetry;
    int formEdges;
    int parmaLoops;
    int parmaVerbosity;
/** \brief write the geombc file during in-memory data transfer
   \ between phasta and chef. */
    int writeGeomBCFiles;
    int ramdisk;
/** \brief the value of criteria for the mesh measure.
   \details this is only used in solver-adaptor (phastaChef) loop.
   \ If the mesh quality is less than this value,
   \ it will trigger the mesh adaptor. */
    double meshqCrtn;
    double elementImbalance;
    double vertexImbalance;
    FILE* (*openfile_read)(Input& in, const char* path);
    RStream* rs;
/** \brief minimum desired mean ratio cubed for simplex elements
   \details a different measure is used for curved elements */
/** \brief the flag for switch between simmetrix mesh and pumi-based mesh.
   \ avoid run incompatible APIs with simmetrix mesh */
    int simmetrixMesh;
/** \brief the max number of iterations in mesh adaptation
   \details this is only used in solver-adaptor (phastaChef) loop */
    int maxAdaptIterations;
    double adaptShrinkLimit;
};

int countNaturalBCs(Input& in);
int countEssentialBCs(Input& in);
int countScalarBCs(Input& in);

}

#endif


