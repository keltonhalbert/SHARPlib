## Add doctest as an interface library
set(DOCTEST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
#add_library(doctest)
include_directories(${CMAKE_SOURCE_DIR}/include ${DOCTEST_INCLUDE_DIR})

## add test executables
add_executable(test_sharp_thermo test_thermo.cpp)
target_link_libraries(test_sharp_thermo SHARP) 

add_executable(test_sharp_interp test_interp.cpp)
target_link_libraries(test_sharp_interp SHARP) 

add_executable(test_sharp_layer test_layer.cpp)
target_link_libraries(test_sharp_layer SHARP) 

add_executable(test_sharp_winds test_winds.cpp)
target_link_libraries(test_sharp_winds SHARP) 

add_executable(test_sharp_parcel test_parcel.cpp)
target_link_libraries(test_sharp_parcel SHARP) 
