function(get_version_from_git)
    find_package(Git QUIET)
    if(NOT Git_FOUND)
        message(WARNING "Git not found")
        return()
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always --match "v[0-9]*"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE GIT_RESULT
    )

    if(NOT GIT_RESULT EQUAL 0)
        message(WARNING "Failed to get git tag")
        return()
    endif()

    # Remove leading 'v'
    string(REGEX REPLACE "^v" "" VERSION_RAW "${GIT_DESCRIBE}")

    # Parse format like: 1.2.3-rc.1-4-gabc1234-dirty
    set(BUILD_METADATA "")
    set(PRE_RELEASE "")
    set(VERSION_STRING "")

    # Check if we're in a dirty state
    if(VERSION_RAW MATCHES "-dirty$")
        set(IS_DIRTY 1)
        string(REGEX REPLACE "-dirty$" "" VERSION_RAW "${VERSION_RAW}")
    else()
        set(IS_DIRTY 0)
    endif()

    # Extract base version, pre-release, build metadata
    if(VERSION_RAW MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)(-([a-zA-Z0-9.]+))?(-([0-9]+)-g([0-9a-f]+))?$")
        set(PROJECT_VERSION_MAJOR "${CMAKE_MATCH_1}")
        set(PROJECT_VERSION_MINOR "${CMAKE_MATCH_2}")
        set(PROJECT_VERSION_PATCH "${CMAKE_MATCH_3}")
        set(PRE_RELEASE "${CMAKE_MATCH_5}")
        set(COMMIT_COUNT "${CMAKE_MATCH_7}")
        set(COMMIT_HASH "${CMAKE_MATCH_8}")

        # Fallback hash if tag is exact
        if(NOT COMMIT_HASH)
            execute_process(
                COMMAND ${GIT_EXECUTABLE} rev-parse --short=7 HEAD
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        endif()

        set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
        if(PRE_RELEASE)
            set(PROJECT_VERSION "${PROJECT_VERSION}-${PRE_RELEASE}")
        endif()

        # Build metadata
        if(COMMIT_COUNT AND COMMIT_HASH)
            set(BUILD_METADATA "${COMMIT_COUNT}.g${COMMIT_HASH}")
        elseif(COMMIT_HASH)
            set(BUILD_METADATA "g${COMMIT_HASH}")
        endif()

        if(IS_DIRTY EQUAL 1)
            if(BUILD_METADATA)
                set(BUILD_METADATA "${BUILD_METADATA}.dirty")
            else()
                set(BUILD_METADATA "dirty")
            endif()
        endif()

        if(BUILD_METADATA)
            set(FULL_VERSION "${PROJECT_VERSION}+${BUILD_METADATA}")
        else()
            set(FULL_VERSION "${PROJECT_VERSION}")
        endif()

        # Export variables
        set(PROJECT_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}" PARENT_SCOPE)
        set(PROJECT_VERSION_MINOR "${PROJECT_VERSION_MINOR}" PARENT_SCOPE)
        set(PROJECT_VERSION_PATCH "${PROJECT_VERSION_PATCH}" PARENT_SCOPE)
        set(PROJECT_VERSION "${PROJECT_VERSION}" PARENT_SCOPE)
        set(FULL_VERSION "${FULL_VERSION}" PARENT_SCOPE)
        set(PRE_RELEASE "${PRE_RELEASE}" PARENT_SCOPE)
        set(COMMIT_COUNT "${COMMIT_COUNT}" PARENT_SCOPE)
        set(COMMIT_HASH "${COMMIT_HASH}" PARENT_SCOPE)
        set(GIT_COMMIT_HASH "g${COMMIT_HASH}" PARENT_SCOPE)
        set(BUILD_METADATA "${BUILD_METADATA}" PARENT_SCOPE)
        set(IS_DIRTY "${IS_DIRTY}" PARENT_SCOPE)
    else()
        message(WARNING "Git tag '${VERSION_RAW}' does not match expected SemVer format.")
    endif()
endfunction()
