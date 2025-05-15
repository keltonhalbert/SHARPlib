## We want to add all of our tests to a custom target and exclude them from "make all" so that
## testing can be built and run optionally. The "make test" target will run the unit tests, but not 
## compile the tests themselves. We create the target SHARPlib_tests so that all unit tests can be
## aggregated under a single build command.
add_custom_target(SHARPlib_tests)

include_directories(${CMAKE_SOURCE_DIR}/include ${DOCTEST_INCLUDE_DIR})

## Create an executable with the given name from the given source, excluding
## it from the "make all" target. Then link it to the base static library,
## and then add a test to CTests
add_executable(test_sharp_thermo EXCLUDE_FROM_ALL test_thermo.cpp)
target_link_libraries(test_sharp_thermo SHARPlib) 
if (NO_QC)
    target_compile_definitions(test_sharp_thermo PRIVATE NO_QC)
endif()
add_test(NAME sharp_thermo COMMAND test_sharp_thermo)

add_executable(test_sharp_interp EXCLUDE_FROM_ALL test_interp.cpp)
target_link_libraries(test_sharp_interp SHARPlib) 
if (NO_QC)
    target_compile_definitions(test_sharp_interp PRIVATE NO_QC)
endif()
add_test(NAME sharp_interp COMMAND test_sharp_interp)

add_executable(test_sharp_layer EXCLUDE_FROM_ALL test_layer.cpp)
target_link_libraries(test_sharp_layer SHARPlib) 
if (NO_QC)
    target_compile_definitions(test_sharp_layer PRIVATE NO_QC)
endif()
add_test(NAME sharp_layer COMMAND test_sharp_layer)

add_executable(test_sharp_winds EXCLUDE_FROM_ALL test_winds.cpp)
target_link_libraries(test_sharp_winds SHARPlib) 
if (NO_QC)
    target_compile_definitions(test_sharp_winds PRIVATE NO_QC)
endif()
add_test(NAME sharp_winds COMMAND test_sharp_winds)

add_executable(test_sharp_parcel EXCLUDE_FROM_ALL test_parcel.cpp)
target_link_libraries(test_sharp_parcel SHARPlib) 
if (NO_QC)
    target_compile_definitions(test_sharp_parcel PRIVATE NO_QC)
endif()
add_test(NAME sharp_parcel COMMAND test_sharp_parcel)

add_executable(test_sharp_algs EXCLUDE_FROM_ALL test_algorithms.cpp)
target_link_libraries(test_sharp_algs SHARPlib) 
if (NO_QC)
    target_compile_definitions(test_sharp_algs PRIVATE NO_QC)
endif()
add_test(NAME sharp_algs COMMAND test_sharp_algs)

## This aggregates all of the executable names into a single 
## target for building and testing
add_dependencies(SHARPlib_tests test_sharp_thermo test_sharp_interp test_sharp_layer test_sharp_winds test_sharp_parcel test_sharp_algs)
