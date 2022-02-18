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
 * @brief get entity imbalance
 * @remark The imbalance of a given entity order (i.e., vtx, edge, face, rgn) 
 * is defined as the maximum count of that entity order on a part, across all parts, 
 * divided by the average entity count across all parts.
 * For example if there are four parts and the parts have 5, 7, 12, and 8
 * vertices, respectively, then the vertex imbalance is 50%;
 * 12 / ((5+7+8+12)/4) = 1.5.
 * @param mesh (InOut) partitioned mesh
 * @param entImb (InOut) entity imbalance [vtx, edge, face, rgn]
 */
void Parma_GetEntImbalance(apf::Mesh* mesh, double (*entImb)[4]);

/**
 * @brief see Parma_GetEntImbalance(...)
 * @remark The weighted imbalance definition replaces the entity count with the 
 * sum of entity weights.  If the weight for all the entities of a given order
 * is one, then the two definitions are equivalent.
 * On a part, if an entity order (vtx, edge, face, rgn) does not have
 * weights set on all its entities then a weight of one will be assigned to each
 * of the entities (of the given order).
 * @param mesh (InOut) partitioned mesh
 * @param weight (In) entity weights used for computing imbalance
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
 * @param maxNumParts (InOut) number of parts with max neighbors
 * @param avg (InOut) average neighbors
 * @param loc (InOut) local neighbors
 */
void Parma_GetNeighborStats(apf::Mesh* m, int& max, int& maxNumParts,
    double& avg, int& loc);

/**
 * @brief write the number of parts with neighbors formed by a small number of shared vtx
 * @param m (In) partitioned mesh
 * @param small (In) report part counts with [1:small] number of shared vertices
 * @param prefix (In) string to prepend to output
 */
void Parma_WriteSmallNeighbors(apf::Mesh* m, int small, const char* prefix);

/**
 * @brief get the smallest number of shared vertices forming a neighbor
 *        ,a 'side', in a part with the maximum number of neigbhors
 * @param m (In) partitioned mesh
 * @return smallest number of shared vertices
 */
int Parma_GetSmallestSideMaxNeighborParts(apf::Mesh* m);

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
 * @brief prints partition stats
 * @remark includes face-disconnected components, number of vertices on
 *         inter-part boundaries, number of vtx-connected neighboring parts,
 *         entity imbalance, and number of empty parts
 * @param m (In) partitioned mesh
 * @param key (In) identifying string to write with stat output
 * @param fine (In) enable per part stat output
 */
void Parma_PrintPtnStats(apf::Mesh* m, std::string key, bool fine=false);

/**
 * @brief prints partition stats using entity weights
 * @remark On a part, if an entity order (vtx, edge, face, rgn) does not have
 * weights set on all its entities then a weight of one will be assigned to each
 * of the entities (of the given order)
 * @param m (In) partitioned mesh
 * @param w (In) tag with entity weights
 * @param key (In) identifying string to write with stat output
 * @param fine (In) enable per part stat output
 */
void Parma_PrintWeightedPtnStats(apf::Mesh* m, apf::MeshTag* w, std::string key, bool fine=false);

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
 * @brief create an APF Balancer for MPAS
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
apf::Balancer* Parma_MakeMPASDiffuser(apf::Mesh* m, int layers, int bridge,
    double stepFactor = 0.1, int verbosity=0);

/**
 * @brief create an edge balancer that is ghost aware
 * @param m (In) partitioned mesh
 * @param stepFactor (In) amount of weight to migrate between parts during
                          diffusion, lower values migrate fewer
                          elements per iteration
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeGhostEdgeDiffuser(apf::Mesh* m,
    double stepFactor = 0.1, int verbosity=0);

/**
 * @brief create an APF Balancer using ghost element aware diffusion for a
 * vertex-based partition
 * @remark Ghosting for a vertex-based partition is asymetric; ghosts
 * are only needed for vertices on the boundary that are owned.
 * @param m (In) partitioned mesh
 * @param layers (In) depth of ghosting
 * @param stepFactor (In) amount of weight to migrate between parts during
                          diffusion, lower values migrate fewer
                          elements per iteration
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m, int layers,
    double stepFactor = 0.1, int verbosity=0);

/** @brief backward compatability
 */
apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m, int layers, int bridge,
    double stepFactor = 0.1, int verbosity=0);

/**
 * @brief write the vertex based partition to file
 * @param m (In) partitioned mesh
 * @param prefix (In) prefix for file names
 */
void Parma_WriteVtxPtn(apf::Mesh* m, const char* prefix);

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

/**
 * @brief reorder the mesh via a breadth first search
 * @remark the returned tag has the reordered vertex order
 * @param m (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf mesh tag
 */
apf::MeshTag* Parma_BfsReorder(apf::Mesh* m, int verbosity=0);

#endif
