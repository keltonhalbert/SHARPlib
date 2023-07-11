
if (${DO_BENCH})
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

	add_subdirectory(../../external/benchmark EXCLUDE_FROM_ALL bench)


	add_custom_target(SHARPlib_bench)

	add_executable(bench_sharp_algorithms EXCLUDE_FROM_ALL bench_algorithms.cpp)
	add_executable(bench_sharp_interp EXCLUDE_FROM_ALL bench_interp.cpp)
	add_executable(bench_sharp_parcel EXCLUDE_FROM_ALL bench_parcel.cpp)

	target_link_libraries(bench_sharp_algorithms SHARPlib benchmark::benchmark)
	target_link_libraries(bench_sharp_interp SHARPlib benchmark::benchmark)
	target_link_libraries(bench_sharp_parcel SHARPlib benchmark::benchmark)

	add_dependencies(SHARPlib_bench bench_sharp_algorithms)
	add_dependencies(SHARPlib_bench bench_sharp_interp)
	add_dependencies(SHARPlib_bench bench_sharp_parcel)

endif()


