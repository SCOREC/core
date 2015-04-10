SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    pcu                 pcu             SS  REQUIRED
    gmi                 gmi             SS  REQUIRED
    gmi_sim             gmi_sim         SS  REQUIRED
    apf                 apf             SS  REQUIRED
    apf_sim             apf_sim         SS  OPTIONAL
    mds                 mds             SS  REQUIRED
    parma               parma           SS  REQUIRED
    apf_zoltan          zoltan          SS  REQUIRED
    ma                  ma              SS  REQUIRED
    spr                 spr             SS  REQUIRED
    apf_stk             stk             SS  REQUIRED
    exe                 test            SS  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS MPI)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
