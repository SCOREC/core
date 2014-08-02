/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef PARMA_H
#define PARMA_H

#include "apf.h"
#include "apfMesh2.h"
#include "apfPartition.h"

/**
 * @brief run partition improvement and heavy part splitting 
 * @param mesh (InOut) partitioned mesh
 * @param weight (In) element weight used for computing imbalance
 * @param maxImb (In) maximum imbalance tolerance
 */
int Parma_Run(apf::Mesh* mesh, apf::MeshTag* weight, const double maxImb);

/**
 * @brief run multi criteria partition improvement
 * @remark see http://redmine.scorec.rpi.edu/projects/parma/wiki for more info
 *
 * @param mesh (InOut) partitioned mesh
 * @param priority (In) entity improvement priority [vtx, edge, face, rgn]
 *                      larger values (>0) have higher priority
 * @param maxImb (In) maximum imbalance tolerance
 * @param verbosity (In) 0: minimal output, >0 increasing amounts
 *                       runtime information
 * @param maxItr (In) maximum diffusion iterations per entity type
 *
 * @return zero on success, non-zero otherwise
 */
int Parma_RunPtnImprovement(apf::Mesh* mesh, int (*priority)[4], 
    const double maxImb=1.05, const int verbosity=0, const int maxItr=20);

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
 * @brief see Parma_RunPtnImprovement(apf::Mesh ... )
 * @param weight (In) element weight used for computing imbalance
 */
int Parma_RunWeightedPtnImprovement(apf::Mesh* mesh, apf::MeshTag* weight, 
    int (*priority)[4], const double maxImb=1.05, const int dbgLvl=0, 
    const int maxItr=20);

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
 * @brief get the maximum and average number of vtx-connected neighboring parts
 * @remark for each part count the number of parts it shares mesh
 *         vertices with
 * @param mesh (In) partitioned mesh
 * @param max (InOut) max neighbors
 * @param avg (InOut) average neighbors
 * @param loc (InOut) local neighbors
 * count
 */
void Parma_GetNeighborStats(apf::Mesh* m, int& max, double& avg, int& loc);

/**
 * @brief get the number of vertices on inter-part boundaries
 * @param mesh (In) partitioned mesh
 * @return number of vertices
 */
long Parma_GetNumBdryVtx(apf::Mesh* m);

/**
 * @brief get the maximum, average and local number of face-disconnected
 * components 
 * @param mesh (In) partitioned mesh
 * @param max (InOut) max disconnected
 * @param avg (InOut) average disconnected
 * @param isMax (InOut) local disconnected 
 */
void Parma_GetDisconnectedStats(apf::Mesh* m, int& max, double& avg, int& loc);

/**
 * @brief prints partition stats
 * @remark includes face-disconnected components, number of vertices on 
 *         inter-part boundaries, number of vtx-connected neighboring parts, 
 *         entity imbalance, and number of empty parts
 * @param mesh (In) partitioned mesh
 * @param key (In) identifying string to write with stat output
 */
void Parma_PrintPtnStats(apf::Mesh* m, std::string key);

/**
 * @brief re-connect disconnected parts
 * @param mesh (In) partitioned mesh
 */
void Parma_ProcessDisconnectedParts(apf::Mesh* m);

/**
 * @brief create an APF Balancer using centroid diffusion
 * @param mesh (In) partitioned mesh
 * @param stepFactor (In) amount of weight to migrate between parts during diffusion, lower values migrate fewer elements per iteration
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeCentroidDiffuser(apf::Mesh* m, double stepFactor = 0.1);

/**
 * @brief create an APF Balancer using ghost element aware diffusion
 * @param mesh (In) partitioned mesh
 * @param layers (In) depth of ghosting
 * @param bridge (In) dimension of entity ghosting depth is based on, typically meshDim-1
 * @param stepFactor (In) amount of weight to migrate between parts during diffusion, lower values migrate fewer elements per iteration
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m, int layers, int bridge, 
    double stepFactor = 0.1, int verbosity=0);

/**
 * @brief create an APF Balancer using heavy part splitting
 * @param mesh (In) partitioned mesh
 * @param verbosity (In) output control, higher values output more
 * @return apf balancer instance
 */
apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, int verbosity=0);

/**
 * @brief create an APF Splitter using recursive inertial bisection
 * @param mesh (In) partitioned mesh
 * @param sync (In) true if all parts will be split, false o.w.
 * @return apf splitter instance
 */
apf::Splitter* Parma_MakeRibSplitter(apf::Mesh* m, bool sync = true);

/**
 * @brief create a mesh tag that weighs elements by their memory consumption
 * @param mesh (In) partitioned mesh
 * @return mesh tag
 */
apf::MeshTag* Parma_WeighByMemory(apf::Mesh* m);

void Parma_ShrinkPartition(apf::Mesh2* m, int factor,
    void (*runAfter)(apf::Mesh2* m));

#endif
