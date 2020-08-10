#ifndef SIM_HELPER_H
#define SIM_HELPER_H


#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>

static bool is_started = false;

void start_sim(const char* logfile = 0);
void stop_sim();
bool is_sim_started();

#endif
