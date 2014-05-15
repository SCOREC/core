SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    pcu                 pcu             SS  REQUIRED
    pumi_util           pumi_util       SS  REQUIRED
    pumi_geom           pumi_geom       SS  REQUIRED
    pumi_mesh           pumi_mesh       SS  REQUIRED
    pumi_geom_meshmodel pumi_geom/meshModel   SS  OPTIONAL
    apf                 apf             SS  OPTIONAL
    apf_pumi            apf/pumi        SS  OPTIONAL
    spr                 apf/spr         SS  OPTIONAL
    apf_stk             apf/stk         SS  OPTIONAL
    parma               parma           SS  OPTIONAL
    meshadapt           meshadapt       SS  OPTIONAL
    ma                  ma              SS  OPTIONAL
)

IF(TPL_ENABLE_Parasolid)
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    ${SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS}
    pumi_geom_parasolid pumi_geom_parasolid   SS  OPTIONAL
    )
endif()
IF(TPL_ENABLE_ACIS)
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    ${SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS}
    pumi_geom_acis      pumi_geom_acis  SS  OPTIONAL
    )
endif()

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS MPI)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

# Turn off the geometry packages by default
set_default(Trilinos_ENABLE_SCORECpumi_geom_parasolid OFF)
set_default(Trilinos_ENABLE_SCORECpumi_geom_acis OFF)

# Turn on meshadapt 2 and parma by default
set_default(Trilinos_ENABLE_SCORECparma ON)
set_default(Trilinos_ENABLE_SCORECma ON)
set_default(Trilinos_ENABLE_SCORECmeshadapt OFF)
