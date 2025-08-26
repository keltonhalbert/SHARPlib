if(NOT STUBS_DIR)
    message(FATAL_ERROR "STUBS_DIR was not defined for fix_nanobind_stubs.cmake")
endif()

file(GLOB SUBMODULE_DIRS LIST_DIRECTORIES true "${STUBS_DIR}/*")

# Filter out any plain .pyi files from the list
list(FILTER SUBMODULE_DIRS EXCLUDE REGEX ".*\\.pyi$")

foreach(SUBMODULE_DIR ${SUBMODULE_DIRS})
    get_filename_component(SUBMODULE_NAME ${SUBMODULE_DIR} NAME)

    set(OLD_PATH "${SUBMODULE_DIR}/__init__.pyi")
    set(NEW_PATH "${STUBS_DIR}/${SUBMODULE_NAME}.pyi")

    message(STATUS "Fixing stub: ${SUBMODULE_NAME}")
    execute_process(COMMAND "${CMAKE_COMMAND}" -E rename "${OLD_PATH}" "${NEW_PATH}")
    execute_process(COMMAND "${CMAKE_COMMAND}" -E remove_directory "${SUBMODULE_DIR}")
endforeach()
