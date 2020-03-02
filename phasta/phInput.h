#ifndef PH_INPUT_H
#define PH_INPUT_H

/** \file phInput.h
    \brief The Chef interface for execution control
    \details The variables defined here should be placed in a file named
            'adapt.inp'.  Each variable should be on its own line and
            followed by a space and then the value assigned to the
            variable. Blank lines and lines starting with '#' are ignored.
            To add a new variable edit this file and phInput.cc .
*/

#include <string>
#include <vector>

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
    /** \brief select the method used to increase the number of parts in the mesh.
        \details partitionMethod can be set to 'graph' to use multi-level
      ParMETIS Part k-way, 'rib' to use SCOREC's recursive inertial bisection,
      and 'zrib' to use Zoltan's recursive inertial bisection. */
    std::string partitionMethod;
    /** \brief select the method used to balance the mesh prior to adaptation.
        \details valid options are 'graph', 'zrib', 'parma', and 'none'.  Selecting
      'parma' balances the elements via a diffusive method and selecting 'none' will
      disable balancing prior to adaptation.  See the partitionMethod documentation
      for a description of the other methods.  Note, the 'LocalPtn' parameter
      does not apply as balancing is a global operation.
     */
    std::string preAdaptBalanceMethod;
    /** \brief select the method used to balance the mesh during adaptation.
        \details valid options are 'graph', 'parma', and 'none'.  See the
      partitionMethod and preAdaptBalanceMethod documentation for a description
      of the methods.  Note, the 'LocalPtn' parameter does not apply as
      balancing is a global operation.
      */
    std::string midAdaptBalanceMethod;
    /** \brief select the method used to balance the mesh after adaptation.
        \details valid options are 'graph', 'zrib', 'parma', 'parma-gap', and 'none'.
      Selecting 'parma-gap' balances the mesh elements and reduces the number of parts
      that share mesh entities with each part (neighbors).  See the partitionMethod
      and preAdaptBalanceMethod documentation for a description of the methods.
      Note, the 'LocalPtn' parameter does not apply as balancing is a global
      operation.
      */
    std::string postAdaptBalanceMethod;
    /** \brief select the method used to balance the mesh prior to pre-processing.
        \details valid options are 'graph', 'zrib', 'parma', 'parma-gap', and 'none'.
      See the partitionMethod, preAdaptBalanceMethod, and postAdaptBalanceMethod
      documentation for a description of the methods.  Note, the 'LocalPtn'
      parameter does not apply as balancing is a global operation.
      */
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
    /** \brief enables the use of local partitioning methods.
        \details when set to '1' each process will run a serial instance of the method
      selected with partitionMethod.  When set to '0' all processes coordinate to run
      a parallel instance of the selected partitioning method.  For part counts over 32Ki,
      the memory requirements of a parallel instance of the 'graph' method typically exceeds
      available memory. */
    int localPtn;
    int recursiveUR;
    int dwalMigration;
    int buildMapping;
    int elementsPerMigration;
    int threaded;
    int initBubbles;
    std::string bubbleFileName;
    int formElementGraph;
    int snap;
    int transferParametric;
    /** \brief enable splitting triangle edges of a prismatic boundary layer stack*/
    int splitAllLayerEdges;
    /** \brief filter out a subset of 3-way periodic matches.
       it also filters out DG interface matches. */
    int filterMatches;
    int axisymmetry;
    int formEdges;
    int parmaLoops;
    int parmaVerbosity;
    /** \brief write the geombc file during in-memory data transfer
       between phasta and chef. */
    int writeGeomBCFiles;
    /** \brief write the restart file during in-memory data transfer
       between phasta and chef. */
    int writeRestartFiles;
    int ramdisk;
    /** \brief the value of criteria for the mesh measure.
        \details this is only used in solver-adaptor (phastaChef) loop.
       If the mesh quality is less than this value,
       it will trigger the mesh adaptor. */
    double meshqCrtn;
    double elementImbalance;
    double vertexImbalance;
    FILE* (*openfile_read)(Input& in, const char* path);
    RStream* rs;
    /** \brief the flag for switch between simmetrix mesh and pumi-based mesh.
       avoid run incompatible APIs with simmetrix mesh */
    int simmetrixMesh;
    /** \brief the max number of iterations in mesh adaptation
        \details this is only used in solver-adaptor (phastaChef) loop */
    int maxAdaptIterations;
    double adaptShrinkLimit;
    /** \brief report the time spent in IO */
    int printIOtime;
    /** \brief flag of writing m2g fields to geomBC files */
    int mesh2geom;
    /** \brief closest distance from zero level set for banded refinement */
    double alphaDist;
    /** \brief absolute isotropic size within [0:alphaDist) of zero level set */
    double alphaSize;
    /** \brief second closest distance from zero level set for banded refinement */
    double betaDist;
    /** \brief absolute isotropic size within [alphaDist:betaDist) of zero level set */
    double betaSize;
    /** \brief furthest distance from zero level set for banded refinement */
    double gammaDist;
    /** \brief absolute isotropic size within [betaDist:gammDist) of zero level set */
    double gammaSize;
    /** \brief number of rigid bodies */
    int nRigidBody;
    /** \brief number of parameters for each rigid body */
    int nRBParam;
    /** \brief parameter data for rigid body */
    std::vector<double> rbParamData;
    /** \brief factor \beta used for mesh smooth/gradation */
    double gradingFactor;
    /** \brief option used for wrapping sim mesh adapter and improver into mover */
    int simCooperation;
    /** \brief flag for writing simmetrix log file */
    int writeSimLog;
    /** \brief maximum CFL number for mesh size */
    double simCFLUpperBound;
    /** \brief minimum desired mesh size for sim adapter */
    double simSizeLowerBound;
    /** \brief maximum desired mesh size for sim adapter */
    double simSizeUpperBound;
    /** \brief number of allowed mesh elements of adapted mesh */
    double simMaxAdaptMeshElements;
};

int countNaturalBCs(Input& in);
int countEssentialBCs(Input& in);
int countScalarBCs(Input& in);

}

#endif


