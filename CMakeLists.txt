# 3.12 for c++20

cmake_minimum_required(VERSION 3.12...3.28)

project(
    rtpghi
    VERSION 1.0
    LANGUAGES C CXX)

# Default build type to Debug.

set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Use ccache if found.

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}") # CMake 3.9+
endif()

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)

add_library(rtpghi SHARED
    src/circularbuf_p.cpp
    src/circularbuf.cpp
    src/pv.cpp
    src/rtdgtreal_p.cpp
    src/rtdgtreal.cpp
    src/rtpghi_heap.cpp
    src/rtpghi_p.cpp
    src/rtpghi.cpp)

target_include_directories(rtpghi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(rtpghi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    PkgConfig::FFTW)

target_link_libraries(rtpghi PRIVATE
    PkgConfig::FFTW)

# Set C++ standard to C++20 (no extensions).

target_compile_features(rtpghi PUBLIC cxx_std_20)
set_target_properties(rtpghi PROPERTIES CXX_EXTENSIONS OFF)

# Enable IPO (Link-Time Optimization) if release build and supported.

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supports_ipo)
    if(supports_ipo)
      set_target_properties(rtpghi PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()
