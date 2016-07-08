/*
 * Copyright 2016 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_FILE_H
#define APF_FILE_H

#include <cstdio>

namespace apf {

class Mesh;

void save_meta(FILE* file, apf::Mesh* mesh);
void restore_meta(FILE* file, apf::Mesh* mesh);

}

#endif
