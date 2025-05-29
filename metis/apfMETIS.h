/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#ifndef APF_METIS_H
#define APF_METIS_H

/**
 * \page metis APF-METIS
 *
 * METIS is a set of serial programs for partitioning graphs. More information
 * is available at https://github.com/KarypisLab/METIS.
 *
 * This module implements a PUMI interface to use METIS as an apf::Splitter or
 * as an apf::Balancer. The Balancer localizes only relevant graph information
 * and tries to minimize migration.
 *
 * The interface is in apfMETIS.h
 */

/**
 * \file apfMETIS.h
 * \brief METIS partitioning for apf::Mesh objects.
 */
namespace apf {

class Mesh;
class Splitter;
class Balancer;

Splitter* makeMETISsplitter(Mesh* mesh);
Balancer* makeMETISbalancer(Mesh* mesh);

} // namespace apf

#endif // APF_METIS_H
