/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRTRILINOSTYPES_H
#define AWRTRILINOSTYPES_H

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>

namespace awr {

typedef int LO;
typedef long GO;
typedef double ST;
typedef Kokkos::DefaultNode::DefaultNodeType NT;
typedef Tpetra::Operator<ST,LO,GO,NT> OP;
typedef Tpetra::Map<LO,GO,NT> MapT;
typedef Tpetra::Vector<ST,LO,GO,NT> VectorT;
typedef Tpetra::CrsMatrix<ST,LO,GO,NT> MatrixT;
typedef Belos::SolverManager<ST,VectorT,OP> SolverT;
typedef Belos::LinearProblem<ST,VectorT,OP> ProblemT;

}

#endif
