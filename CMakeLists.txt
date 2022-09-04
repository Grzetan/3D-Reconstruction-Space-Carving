cmake_minimum_required(VERSION 3.9)
project(Polyhedron)

#set(CMAKE_CXX_COMPILER_ID "MSVC")


if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)

set(CMAKE_VERBOSE_MAKEFILE ON)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_FLAGS "-std=c++17")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    if (MSVC_VERSION GREATER_EQUAL "1900")
        include(CheckCXXCompilerFlag)
        CHECK_CXX_COMPILER_FLAG("/std:c++latest" _cpp_latest_flag_supported)
        if (_cpp_latest_flag_supported)
            add_compile_options("/std:c++latest")
        endif ()
    endif ()
endif ()


set(PATH ./src)

set(SOURCE_FILES
        ${PATH}/happly.h
        ${PATH}/params.h
        ${PATH}/params.cpp
        ${PATH}/types.h
        ${PATH}/utils.h
        ${PATH}/utils.cpp
        )

add_executable(main ${PATH}/main.cpp ${SOURCE_FILES})


