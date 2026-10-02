## We want to add all of our tests to a custom target and exclude them from "make all" so that
## testing can be built and run optionally. The "make test" target will run the unit tests, but not 
## compile the tests themselves. We create the target SHARPlib_tests so that all unit tests can be
## aggregated under a single build command.
add_custom_target(SHARPlib_tests)

include_directories(${CMAKE_SOURCE_DIR}/include ${DOCTEST_INCLUDE_DIR})

## Fails to compile if NO_QC stops propagating from the SHARPlib target
add_executable(test_sharp_no_qc_propagation EXCLUDE_FROM_ALL test_no_qc_propagation.cpp)
target_link_libraries(test_sharp_no_qc_propagation SHARPlib)
target_compile_definitions(test_sharp_no_qc_propagation PRIVATE SHARPLIB_EXPECT_NO_QC=$<BOOL:${NO_QC}>)
add_test(NAME sharp_no_qc_propagation COMMAND test_sharp_no_qc_propagation)

## Create an executable with the given name from the given source, excluding
## it from the "make all" target. Then link it to the base static library,
## and then add a test to CTests
foreach(suite thermo interp layer winds parcel params qc)
    add_executable(test_sharp_${suite} EXCLUDE_FROM_ALL test_${suite}.cpp)
    target_link_libraries(test_sharp_${suite} SHARPlib)
    add_test(NAME sharp_${suite} COMMAND test_sharp_${suite})
    list(APPEND SHARP_TEST_TARGETS test_sharp_${suite})
endforeach()

add_executable(test_sharp_algs EXCLUDE_FROM_ALL test_algorithms.cpp)
target_link_libraries(test_sharp_algs SHARPlib) 
add_test(NAME sharp_algs COMMAND test_sharp_algs)

## This aggregates all of the executable names into a single 
## target for building and testing
add_dependencies(SHARPlib_tests ${SHARP_TEST_TARGETS} test_sharp_algs test_sharp_no_qc_propagation)
