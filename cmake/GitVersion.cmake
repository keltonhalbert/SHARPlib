# In cmake/GitVersion.cmake or wherever you define get_version_from_git
function(get_project_version)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
        # --- This is your original logic for developer builds ---
        find_package(Git QUIET)
        if(NOT Git_FOUND)
            message(FATAL_ERROR "Git not found, but .git directory exists.")
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always --match "v[0-9]*"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_DESCRIBE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE GIT_RESULT
        )
        if(NOT GIT_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to get git tag for a developer build.")
        endif()

        # We can reuse part of your parsing logic.
        # This simplified version just gets the version for the header.
        string(REGEX REPLACE "^v" "" FULL_VERSION "${GIT_DESCRIBE}")
        string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" PROJECT_VERSION "${FULL_VERSION}")
        string(REGEX MATCH "^([0-9]+)" PROJECT_VERSION_MAJOR "${PROJECT_VERSION}")
        string(REGEX MATCH "\\.([0-9]+)" PROJECT_VERSION_MINOR "${PROJECT_VERSION}")
        set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_1}) # Extract from ".<minor>"
        string(REGEX MATCH "\\.([0-9]+)$" PROJECT_VERSION_PATCH "${PROJECT_VERSION}")
        set(PROJECT_VERSION_PATCH ${CMAKE_MATCH_1}) # Extract from ".<patch>"

        # You can add more detailed parsing here if needed for the C++ header
        set(PRE_RELEASE "" PARENT_SCOPE)
        set(BUILD_METADATA "" PARENT_SCOPE)
        set(IS_DIRTY 0 PARENT_SCOPE)
        if(FULL_VERSION MATCHES "-dirty")
            set(IS_DIRTY 1 PARENT_SCOPE)
        endif()

    else()
        # --- This is the new logic for sdist builds ---
        message(STATUS "No .git directory found. Reading version from _version.py for sdist build.")
        set(VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/nanobind/_version.py")
        if(NOT EXISTS "${VERSION_FILE}")
            message(FATAL_ERROR "sdist build failed: ${VERSION_FILE} not found.")
        endif()

        file(READ ${VERSION_FILE} VERSION_PY_CONTENTS)
        
        # CORRECTED REGEX: Looks for "version = '...'" or "version = \"...\""
        string(REGEX MATCH "version = ['\"]([^'\"]+)['\"]" _ "${VERSION_PY_CONTENTS}")

        if(NOT CMAKE_MATCH_1)
            message(FATAL_ERROR "Could not parse version from ${VERSION_FILE}")
        endif()
        set(FULL_VERSION "${CMAKE_MATCH_1}")

        # Parse the version string for the header file
        string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" PROJECT_VERSION "${FULL_VERSION}")
        set(PROJECT_VERSION_MAJOR "${CMAKE_MATCH_1}")
        set(PROJECT_VERSION_MINOR "${CMAKE_MATCH_2}")
        set(PROJECT_VERSION_PATCH "${CMAKE_MATCH_3}")

        # In an sdist, there's no extra metadata
        set(PRE_RELEASE "" PARENT_SCOPE)
        set(BUILD_METADATA "" PARENT_SCOPE)
        set(IS_DIRTY 0 PARENT_SCOPE)
    endif()

    # Export variables to the parent scope for version.h and project() command
    set(PROJECT_VERSION "${PROJECT_VERSION}" PARENT_SCOPE)
    set(FULL_VERSION "${FULL_VERSION}" PARENT_SCOPE)
    set(PROJECT_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}" PARENT_SCOPE)
    set(PROJECT_VERSION_MINOR "${PROJECT_VERSION_MINOR}" PARENT_SCOPE)
    set(PROJECT_VERSION_PATCH "${PROJECT_VERSION_PATCH}" PARENT_SCOPE)

endfunction()
