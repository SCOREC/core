/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apfPartition.h>

namespace apf {

namespace metis {

class MetisSplitter : public Splitter {
public:
  MetisSplitter(Mesh* mesh): mesh_(mesh) {}
  ~MetisSplitter() {}
  Migration* split(MeshTag* weights, double tolerance, int multiple);
private:
  Mesh* mesh_;
};

} // namespace metis

} // namespace apf

