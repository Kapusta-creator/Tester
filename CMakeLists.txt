cmake_minimum_required(VERSION 3.14)

project(tester)

include_directories(prj.lab)
include(FetchContent)

FetchContent_Declare(matplotplusplus GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus GIT_TAG origin/master)

FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
add_subdirectory(prj.lab)
