/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_MIXED_NUMBERING_H
#define APF_MIXED_NUMBERING_H

/** \file apfMixedNumbering.h
  * \brief Global numbering interface for mixed fields.
  * \details This is separate from apfNumbering.h because it pulls
  * in std::vector. */

#include "apf.h"
#include <vector>

namespace apf {

typedef NumberingOf<long> GlobalNumbering;

/** \brief Count the total numbered degrees of freedom.
  * \param n The input local numberings.
  * \returns The number of on-part numbered degrees of freedom */
int countDOFs(std::vector<Numbering*> const& n);

/** \brief Get the element numbers for multiple numberings.
  * \param n The input numberings.
  * \param e The mesh entity for which to get element numbers.
  * \param numbers the output element numbers. */
void getElementNumbers(
    std::vector<Numbering*> const& n,
    MeshEntity* e,
    std::vector<int>& numbers);

/** \brief Get the element numbers for multiple global numberings.
  * \param n The input global numberings.
  * \param e The mesh entity for which to get element numbers.
  * \param numbers The output element numbers. */
void getElementNumbers(
    std::vector<GlobalNumbering*> const& n,
    MeshEntity* e,
    std::vector<long>& numbers);

/** \brief Number the owned nodes of multiple fields.
  * \param fields The input fields to be numbered.
  * \param owned The output local owned numberings.
  * \returns The number of on-part owned nodes across all fields. */
int numberOwned(
    std::vector<Field*> const& fields,
    std::vector<Numbering*>& owned);

/** \brief Number the ghost (overlapped/shared) nodes of multiple fields.
  * \param fields The input fields to be numbered.
  * \param ghost The output local ghost numberings.
  * \returns The number of on-part ghost nodes across all fields. */
int numberGhost(
    std::vector<Field*> const& fields,
    std::vector<Numbering*>& ghost);

/** \brief Globalize a mixed numbering scheme.
  * \param owned The local owned on-part numberings.
  * \param global The output global numberings. */
void makeGlobal(
    std::vector<Numbering*>& owned,
    std::vector<GlobalNumbering*>& global);

}

#endif
