/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apfPartition.h>

namespace apf {

namespace metis {

class MetisBalancer : public Balancer {
public:
  MetisBalancer(Mesh* mesh) : mesh_(mesh) {}
  ~MetisBalancer() {}
  void balance(MeshTag* weights, double tolerance);
private:
  Mesh* mesh_;
};

} // namespace metis

} // namespace apf

