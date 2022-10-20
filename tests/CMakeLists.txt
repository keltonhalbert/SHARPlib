## Add doctest as an interface library
set(DOCTEST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
#add_library(doctest)
include_directories(${CMAKE_SOURCE_DIR}/include ${DOCTEST_INCLUDE_DIR})

## add test executable
add_executable(test_sharp_thermo test_thermo.cpp)
target_link_libraries(test_sharp_thermo SHARP) 
