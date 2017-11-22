/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_STATS_H
#define MA_STATS_H

/* #include <apfMesh2.h> */
/* #include <apfShape.h> */
/* #include <apfNumbering.h> */
#include "maAdapt.h"
#include "maShape.h"
#include <pcu_util.h>

#include <vector>
#include <math.h>

namespace ma {

void getLinearQualitiesInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &linearQualities);

void getEdgeLengthsInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths);

void getLinearQualitiesInPhysicalSpace(ma::Mesh* m,
    std::vector<double> &linearQualities);

void getEdgeLengthsInPhysicalSpace(ma::Mesh* m,
    std::vector<double> &edgeLengths);

void getStatsInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities);

void getStatsInPhysicalSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities);

void stats(ma::Input* in,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    bool inMetric);

}
#endif
