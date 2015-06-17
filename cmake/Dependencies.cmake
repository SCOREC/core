SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    pcu                 pcu             PT  REQUIRED
    gmi                 gmi             PT  REQUIRED
    apf                 apf             PT  REQUIRED
    mds                 mds             PT  REQUIRED
    parma               parma           PT  REQUIRED
    apf_zoltan          zoltan          PT  OPTIONAL
    ma                  ma              PT  OPTIONAL
    spr                 spr             PT  REQUIRED
    apf_stk             stk             PT  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)

SET(LIB_OPTIONAL_DEP_TPLS MPI)

SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
