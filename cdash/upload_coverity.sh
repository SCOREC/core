cd /fasttmp/seol/scorec/cdash/build/core
curl --form token=faMZVmxTjByhNoJyb_4wDw \
  --form email=seols@rpi.edu \
  --form file=@$PWD/pumi.tgz \
  --form version="1.0.1" \
  --form description="Automated" \
  https://scan.coverity.com/builds?project=SCOREC%2Fcore
