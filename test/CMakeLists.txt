find_package(gmi PATHS ${CMAKE_BINARY_DIR})
set(PUMI_INCLUDES ${GMI_INCLUDE_DIRS})
set(PUMI_LIBS ${GMI_LIBS})

include_directories(${PUMI_INCLUDES})

add_executable(gmi_test EXCLUDE_FROM_ALL gmi_test.c)

set(tests gmi_test)

foreach(test ${tests})
  target_link_libraries(${test} ${PUMI_LIBS})
endforeach()
