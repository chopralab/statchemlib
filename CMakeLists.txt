project(statchemlib)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

set(STCH_BUILD_STATIC_EXECUTABLE OFF CACHE BOOL "Build a version of candock without dynamic linking")

if(${STCH_BUILD_STATIC_EXECUTABLE})
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif()
endif()

if (WIN32)

    set(Boost_USE_STATIC_LIBS        ON)
    set(Boost_USE_MULTITHREADED      ON)    

    add_definitions( /D_USE_MATH_DEFINES /DWIN32_LEAN_AND_MEAN /DWINVER=0x0A00 /D_WIN32_WINNT=0x0A00)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC -Wall -Wextra -pedantic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -fPIC -Wall -Wextra -pedantic")
endif()

find_package(Boost    REQUIRED  COMPONENTS filesystem system program_options date_time)
find_package(LibLZMA  REQUIRED)
find_package(Git      REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

find_package(OpenMM REQUIRED COMPONENTS openmm)
find_package(GSL REQUIRED)

include_directories(${OPENMM_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${LIBLZMA_INCLUDE_DIRS})

install(DIRECTORY include/ DESTINATION include)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
install(DIRECTORY data/ DESTINATION data)

# This order is important as CMAKE needs to find certain files
add_subdirectory(lib)
add_subdirectory(src)

option(ENABLE_TESTING "Enable testing for CANDIY projects" OFF)

if (${ENABLE_TESTING})
    enable_testing()
    add_subdirectory(tests)
endif()
