# 3.12 for c++20

cmake_minimum_required(VERSION 3.12...3.28)

project(
    rtpghi
    VERSION 1.0
    LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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
pkg_search_module(sndfile REQUIRED sndfile IMPORTED_TARGET)

add_library(rtpghi STATIC
    include/firwin.h
    include/pv.h
    include/rtdgtreal.h
    include/rtdgtrealproc.h
    include/rtpghi.h
    src/arrayutils.cpp
    src/arrayutils.h
    src/circularbuf.cpp
    src/circularbuf.h
    src/firwin.cpp
    src/gabdual_painless.cpp
    src/gabdual_painless.h
    src/pv_p.cpp
    src/pv_p.h
    src/pv.cpp
    src/rtdgtreal_p.cpp
    src/rtdgtreal_p.h
    src/rtdgtreal.cpp
    src/rtdgtrealproc_p.cpp
    src/rtdgtrealproc_p.h
    src/rtdgtrealproc.cpp
    src/rtpghi_heap.cpp
    src/rtpghi_heap.h
    src/rtpghi_p.cpp
    src/rtpghi_p.h
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

# Test.

add_executable(test tests/main.cpp)

target_link_libraries(test PRIVATE
    rtpghi
    PkgConfig::sndfile)

target_include_directories(test PRIVATE
    PkgConfig::sndfile)

# Enable Sanitizers if debug build.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(rtpghi PUBLIC -fsanitize=address)
    target_link_options(rtpghi PUBLIC -fsanitize=address)
endif()