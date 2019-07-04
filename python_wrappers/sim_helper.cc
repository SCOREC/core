#include<sim_helper.h>

void start_sim(const char* logfile)
{
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  if (logfile)
    Sim_logOn(logfile);
}

void stop_sim()
{
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
}

