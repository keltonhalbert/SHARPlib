 
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

    set(VERSION_RAW "${GIT_DESCRIBE}")

    # 1. Handle dirty state first and remove the suffix
    if(VERSION_RAW MATCHES "-dirty$")
        set(IS_DIRTY 1)
        string(REGEX REPLACE "-dirty$" "" VERSION_RAW "${VERSION_RAW}")
    else()
        set(IS_DIRTY 0)
    endif()

    # 2. Extract commit info (-<count>-g<hash>) if it exists
    if(VERSION_RAW MATCHES "-([0-9]+)-g([0-9a-f]+)$")
        set(COMMIT_COUNT "${CMAKE_MATCH_1}")
        set(COMMIT_HASH "${CMAKE_MATCH_2}")
        # Remove the commit info from the string to isolate the tag part
        string(REGEX REPLACE "-[0-9]+-g[0-9a-f]+$" "" VERSION_RAW "${VERSION_RAW}")
    else()
        set(COMMIT_COUNT 0)
        # If no count/hash, we're on a tag or it's just a hash. Get the hash directly.
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short=7 HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()

    # 3. Parse the remaining tag part (v1.2.3 or v1.2.3-rc1)
    string(REGEX REPLACE "^v" "" VERSION_RAW "${VERSION_RAW}")
    if(VERSION_RAW MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)(-(.*))?$")
        set(PROJECT_VERSION_MAJOR "${CMAKE_MATCH_1}")
        set(PROJECT_VERSION_MINOR "${CMAKE_MATCH_2}")
        set(PROJECT_VERSION_PATCH "${CMAKE_MATCH_3}")
        set(PRE_RELEASE_TAG_PART "${CMAKE_MATCH_5}") # Pre-release from the tag itself
    else()
        # Fallback if no M.m.p tag was found (describe returned only a hash)
        message(WARNING "Could not parse a SemVer tag from '${GIT_DESCRIBE}'. Defaulting to 0.0.1.")
        set(PROJECT_VERSION_MAJOR 0)
        set(PROJECT_VERSION_MINOR 0)
        set(PROJECT_VERSION_PATCH 1)
        set(PRE_RELEASE_TAG_PART "dev")
    endif()

    # 4. Apply the "next version" logic
    if(COMMIT_COUNT GREATER 0)
        # If we are on a commit after a tag, this is a development build.
        set(PRE_RELEASE "dev${COMMIT_COUNT}")

        # ONLY increment the patch number if the base tag was a FINAL release
        # (i.e., it did not have a pre-release part like "-alpha.0").
        if(NOT PRE_RELEASE_TAG_PART)
            math(EXPR PROJECT_VERSION_PATCH "${PROJECT_VERSION_PATCH} + 1")
        endif()
    else()
        # If we are exactly on a tag, use its pre-release component.
        set(PRE_RELEASE "${PRE_RELEASE_TAG_PART}")
    endif()

    # 5. Assemble the final version strings
    set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
    if(PRE_RELEASE)
        set(PROJECT_VERSION "${PROJECT_VERSION}-${PRE_RELEASE}")
    endif()

    # Build metadata
    if(COMMIT_HASH)
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
    set(COMMIT_HASH "${COMMIT_HASH}" PARENT_SCOPE)
    set(GIT_COMMIT_HASH "g${COMMIT_HASH}" PARENT_SCOPE)
    set(BUILD_METADATA "${BUILD_METADATA}" PARENT_SCOPE)
    set(IS_DIRTY "${IS_DIRTY}" PARENT_SCOPE)
endfunction()
