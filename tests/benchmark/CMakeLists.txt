
## Here we are doing some checking for our git submodules and
## downloading/initializing them accordingly. If the submodules
## have not already been recursively downloaded, this will do it.
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
# Update submodules as needed
	option(GIT_SUBMODULE "Check submodules during build" ON)
	if(GIT_SUBMODULE)
		message(STATUS "Submodule update")
		execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
		if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        	message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/benchmark/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()



## Add the benchmark library as a dependecy and 
## set some build options
set(BENCHMARK_ENABLE_TESTING NO)
set(CMAKE_BUILD_TYPE "Release")
add_subdirectory(../../external/benchmark bench EXCLUDE_FROM_ALL)
message(STATUS ${CMAKE_BUILD_TYPE})

add_executable(bench_sharp_algorithms bench_algorithms.cpp)
target_link_libraries(bench_sharp_algorithms benchmark::benchmark)


