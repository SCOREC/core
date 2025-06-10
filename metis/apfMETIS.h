/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#ifndef APF_METIS_H
#define APF_METIS_H

#define APF_METIS_MAXRANKS 256

/**
 * \file apfMETIS.h
 * \brief METIS partitioning for apf::Mesh objects.
 */
namespace apf {

class Mesh;
class Splitter;
class Balancer;

/**
 * \addtogroup Partitioning
 * \{
 */

/**
 * \brief Make an apf::Splitter that calls METIS for the underlying algorithm
 *
 * \param mesh the apf::Mesh to split
 * \return an apf::Splitter object
 */
Splitter* makeMETISsplitter(Mesh* mesh);
/**
 * \brief Make an apf::Balancer that calls METIS for the underlying algorithm
 *
 * The necessary graph information is localized to rank 0 from mesh topology,
 * METIS is run serially, then the 'migration plan' is sent to the remaining
 * N-1 processes.
 *
 * \param mesh the apf::Mesh to balance
 * \return an apf::Balancer object
 */
Balancer* makeMETISbalancer(Mesh* mesh);
/**
 * \brief Query whether METIS is supported.
 *
 * Support for METIS may or may not be included during compilation. This
 * function allows users to check whether the APF-METIS routines will succeed.
 * If this function returns false, all APF-METIS functions will fail.
 *
 * \return true if APF-METIS is compiled with METIS support
 */
bool hasMETIS();

/**\}*/

} // namespace apf

#endif // APF_METIS_H
