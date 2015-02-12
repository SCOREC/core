/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef PARMA_H
#define PARMA_H

/** \file parma.h
    \brief The ParMA tools interface */

#include "apf.h"
#include "apfPartition.h"

/**
 * @brief run ghost partition improvement
 * @param mesh (InOut) partitioned mesh
 * @param weight (In) element weight used for computing imbalance
 * @param maxImb (In) maximum imbalance tolerance
 * @param numlayers (In) number of ghost layers
 * @param bridgeDim (In) dimension of bridge entity
 * @param verbosity (In) 0: minimal output, >0 increasing amounts
 *                       runtime information
 * @return zero on success, non-zero otherwise
 */
int Parma_RunGhostPtnImprovement(apf::Mesh* mesh, apf::MeshTag* weight,
    const double maxImb, const int numlayers, const int bridgeDim,
    const int verbosity=0);

/**
 * @brief run selp to create a partition with kN parts where k is the
 *        partition factor, and N is the number of parts in the input mesh
 * @param mesh (InOut) partitioned mesh with N parts
 * @param weight (In) element weight used for computing imbalance
 * @param factor (In) partition factor > 1
 * @param verbosity (In) 0: minimal output, >0 increasing amounts
 * @return zero on success, non-zero otherwise
 */
int Parma_RunSelp(apf::Mesh* mesh, apf::MeshTag* weight, const int factor,
    const int verbosity=0);

/**
 * @brief get entity imbalance
 * @param mesh (InOut) partitioned mesh
 * @param entImb (InOut) entity imbalance [vtx, edge, face, rgn]
 */
void Parma_GetEntImbalance(apf::Mesh* mesh, double (*entImb)[4]);

/**
 * @brief see Parma_GetEntImbalance(...)
 * @param mesh (InOut) partitioned mesh
 * @param weight (In) element weight used for computing imbalance
 * @param entImb (InOut) entity imbalance [vtx, edge, face, rgn]
 */
void Parma_GetWeightedEntImbalance(apf::Mesh* mesh, apf::MeshTag* weight,
    double (*entImb)[4]);

/**
 * @brief see Parma_GetEntImbalance(...)
 * @param mesh (InOut) partitioned mesh
 * @param weight (In) element weight used for computing imbalance
 * @param dim (In) entity dimension [vtx|edge|face|rgn]
 * @return entity imbalance
 */
double Parma_GetWeightedEntImbalance(apf::Mesh* mesh, apf::MeshTag* weight,
    int dim);


/**
 * @brief get the maximum and average number of vtx-connected neighboring parts
 * @remark for each part count the number of parts it shares mesh
 *         vertices with
 * @param m (In) partitioned mesh
 * @param max (InOut) max neighbors
 * @param avg (InOut) average neighbors
 * @param loc (InOut) local neighbors
 */
void Parma_GetNeighborStats(apf::Mesh* m, int& max, double& avg, int& loc);

/**
 * @brief get the number of owned vertices on inter-part boundaries
 * @param m (In) partitioned mesh
 * @param loc (InOut) local number of vertices
 * @param tot (InOut) total number of vertices
 * @param min (InOut) min number of vertices on a single part
 * @param max (InOut) max number of vertices on a single part
 * @param avg (InOut) average number of vertices per part
 */
void Parma_GetOwnedBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg);

/**
 * @brief get the number of shared vertices on inter-part boundaries
 * @param m (In) partitioned mesh
 * @param loc (InOut) local number of vertices
 * @param tot (InOut) total number of vertices
 * @param min (InOut) min number of vertices on a single part
 * @param max (InOut) max number of vertices on a single part
 * @param avg (InOut) average number of vertices per part
 */
void Parma_GetSharedBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg);

/**
 * @brief get the number of vertices classified on the model boundary
 * @param m (In) partitioned mesh
 * @param loc (InOut) local number of vertices
 * @param tot (InOut) total number of vertices
 * @param min (InOut) min number of vertices on a single part
 * @param max (InOut) max number of vertices on a single part
 * @param avg (InOut) average number of vertices per part
 */
void Parma_GetMdlBdryVtxStats(apf::Mesh* m, int& loc, long& tot, int& min,
    int& max, double& avg);

/**
 * @brief get the maximum, average and local number of face-disconnected
 * components
 * @param m (In) partitioned mesh
 * @param max (InOut) max disconnected
 * @param avg (InOut) average disconnected
 * @param loc (InOut) local disconnected
 */
void Parma_GetDisconnectedStats(apf::Mesh* m, int& max, double& avg, int& loc);

/**
 * @brief get the maximum, average and local number of entities of the 
 *        specified order/dim
 * @param m (In) partitioned mesh
 * @param dim (In) entity order/dimension of interest
 * @param tot (InOut) total ents
 * @param min (InOut) min ents
 * @param max (InOut) max ents
 * @param avg (InOut) average ents
 * @param loc (InOut) local ents
 */
void Parma_GetEntStats(apf::Mesh* m, int dim, long& tot, int& min, int& max, 
    double& avg, int& loc);

/**
 * @brief prints partition stats
 * @remark includes face-disconnected components, number of vertices on
 *         inter-part boundaries, number of vtx-connected neighboring parts,
 *         entity imbalance, and number of empty parts
 * @param m (In) partitioned mesh
 * @param key (In) identifying string to write with stat output
 */
void Parma_PrintPtnStats(apf::Mesh* m, std::string key, bool fine=false);

/**
 * @brief re-connect disconnected parts
 * @param m (In) partitioned mesh
 */
void Parma_ProcessDisconnectedParts(apf::Mesh* m);

/**
 * @brief create an APF Balancer using centroid diffusion
 * @param m (In) partitioned mesh
 * @param stepFactor (In) amount of weight to migrate between parts
                          during diffusion, lower values migrate
                          fewer elements per iteration
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeCentroidDiffuser(apf::Mesh* m, double stepFactor = 0.1,
    int verbose=0);

/**
 * @brief create an APF Balancer to optimize part shape
 * @param m (In) partitioned mesh
 * @param stepFactor (In) amount of weight to migrate between parts during
               diffusion, lower values migrate fewer elements per iteration
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m, double stepFactor = 0.1,
    int verbose=0);

/**
 * @brief create an APF Balancer to weld small part boundaries together
 * @param m (In) partitioned mesh
 * @param stepFactor (In) amount of weight to migrate between parts during
               diffusion, lower values migrate fewer elements per iteration
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeWelder(apf::Mesh* m, double stepFactor = 0.1,
    int verbose=0);

/**
 * @brief create an APF Balancer using ghost element aware diffusion
 * @param m (In) partitioned mesh
 * @param layers (In) depth of ghosting
 * @param bridge (In) dimension of entity ghosting depth is based on,
                      typically meshDim-1
 * @param stepFactor (In) amount of weight to migrate between parts during
                          diffusion, lower values migrate fewer
                          elements per iteration
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m, int layers, int bridge,
    double stepFactor = 0.1, int verbosity=0);

/**
 * @brief create an APF Balancer using heavy part splitting
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, int verbosity=0);

/**
 * @brief create an APF Balancer targeting vertex imbalance
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeVtxBalancer(apf::Mesh* m, double stepFactor=0.1,
    int verbosity=0);

/**
 * @brief create an APF Balancer targeting element imbalance
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeElmBalancer(apf::Mesh* m, double stepFactor=0.1,
    int verbosity=0);

/**
 * @brief create an APF Balancer targeting element imbalance while preserving
 *        vertex balance
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeElmLtVtxBalancer(apf::Mesh* m, double maxVtx,
    double stepFactor=0.1, int verbosity=0);

/**
 * @brief create an APF Balancer targeting vertex, edge, and elm imbalance
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeVtxEdgeElmBalancer(apf::Mesh* m,
    double stepFactor=0.1, int verbosity=0);

/**
 * @brief create an APF Balancer targeting vertex, and elm imbalance
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeVtxElmBalancer(apf::Mesh* m,
    double stepFactor=0.1, int verbosity=0);

/**
 * @brief create an APF Balancer targeting edge imbalance
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeEdgeBalancer(apf::Mesh* m, double stepFactor=0.1,
    int verbosity=0);

/**
 * @brief create an APF Splitter using recursive inertial bisection
 * @param m (In) partitioned mesh
 * @param sync (In) true if all parts will be split, false o.w.
 * @return apf splitter instance
 */
apf::Splitter* Parma_MakeRibSplitter(apf::Mesh* m, bool sync = true);

/**
 * @brief create a mesh tag that weighs elements by their memory consumption
 * @param m (In) partitioned mesh
 * @return mesh tag
 */
apf::MeshTag* Parma_WeighByMemory(apf::Mesh* m);

/**
 * @brief User-defined code to run on process sub-groups.
 */
struct Parma_GroupCode
{
  /**
   * @brief Called withing sub-groups.
   * @details Within a group, all PCU functions behave
   *          as though only the group processes exist.
   *          PCU_Comm_Peers and PCU_Comm_Self reflect
   *          the number of processes in the group and
   *          the process id within the group, and all
   *          collective operations are confined inside
   *          the group.
   * @param group the group id number, starting from zero
   */
  virtual void run(int group) = 0;
};

/**
 * @brief Shrink the mesh into N/factor processes.
 * @details This function will take N=PCU_Comm_Peers() and generate
 *          (factor) subgroups, each of (N/factor) processes.
 *          It then migrates the mesh onto group 0 and then calls
 *          the user's group code on all groups.
 *          After the user's code completes, the mesh is repartitioned
 *          to all N processes and this function returns.
 *          groups are organized such that contiguous ranges of (factor)
 *          parts are combined into one and then that one is split back
 *          out into (factor) contiguous part ids again.
 */
void Parma_ShrinkPartition(apf::Mesh2* m, int factor, Parma_GroupCode& toRun);

/**
 * @brief Split the processes into groups of (factor).
 * @details This function groups contiguous ranges of (factor)
 *          processes into (PCU_Comm_Peers()/factor) total groups,
 *          then calls the user's group code.
 *          After the user's code completes, execution returns
 *          to the global communicator.
 *          If a non-zero mesh pointer is given, then
 *          apf::remapPartition is used to maintain the mesh structure
 *          during these transitions.
 */
void Parma_SplitPartition(apf::Mesh2* m, int factor, Parma_GroupCode& toRun);

/**
 * @brief Compute maximal independent set numbering
 * @remark This function will compute the maximal independent set numbering
 *         for the partition such that no two part neighbors that share
 *         dimension d mesh entities will be assigned the same number.
 * @param m (In) partitioned mesh
 * @param d (In) adjacency dimension
 */
int Parma_MisNumbering(apf::Mesh* m, int d);


#endif
