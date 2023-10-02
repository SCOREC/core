/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ZOLTAN_CALLBACKS_H
#define APF_ZOLTAN_CALLBACKS_H

#include <zoltan.h>

namespace apf {

class ZoltanMesh;

class ZoltanData {
  public:
    ZoltanData(ZoltanMesh* zb_);
    ~ZoltanData();
    void run();
    int getNumExported() {return num_exported;};
    void getExport(int ind, int* localId, int* export_part);
  private:
    ZoltanData();
    void setup();
    void ptn();
    ZoltanMesh* zb;
    Zoltan_Struct* ztn;
    ZOLTAN_ID_PTR import_gids;
    ZOLTAN_ID_PTR import_lids; /* Pointers to nodes imported */
    int *import_procs; /* Proc IDs of procs owning nodes to be imported.*/
    ZOLTAN_ID_PTR export_gids; /* Global ids of nodes exported */
    ZOLTAN_ID_PTR export_lids; /* Pointers to nodes exported */
    int *export_procs;
    int num_imported; /* Number of nodes to be imported. */
    int num_exported; /* Number of nodes to be exported. */
    int *import_to_part;
    int *export_to_part;
    int dbgLvl;
    int changes; 
    int lidSz;
    int gidSz;
};

}

#endif
