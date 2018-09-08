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

find_package(Boost    REQUIRED  COMPONENTS filesystem system program_options)
find_package(LibLZMA  REQUIRED)
find_package(Git      REQUIRED)

if (WIN32)

    set(Boost_USE_STATIC_LIBS        OFF)
    set(Boost_USE_MULTITHREADED      ON)    

	if (${STCH_BUILD_STATIC_EXECUTABLE})
	    set(Boost_USE_STATIC_LIBS        ON)
	    set(Boost_USE_STATIC_RUNTIME     ON)
	    set(BOOST_ALL_DYN_LINK           OFF)   # force dynamic linking for all libraries
	else()
	    set(Boost_USE_STATIC_LIBS        OFF)
	    set(Boost_USE_STATIC_RUNTIME     OFF)
	    set(BOOST_ALL_DYN_LINK           ON)   # force dynamic linking for all libraries
	endif()

    add_definitions( /D_USE_MATH_DEFINES /DBOOST_ALL_NO_LIB /DBOOST_ALL_DYN_LINK /DWIN32_LEAN_AND_MEAN /DWINVER=0x0A00 /D_WIN32_WINNT=0x0A00)

else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC -Wall -Wextra -pedantic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -fPIC -Wall -Wextra -pedantic")
endif()

include_directories(${Boost_INCLUDE_DIRS})

find_package(OpenMM REQUIRED COMPONENTS openmm)
find_package(GSL REQUIRED)

include_directories(${OPENMM_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${LIBLZMA_INCLUDE_DIRS})

install(DIRECTORY include/ DESTINATION include)
include_directories(BEFORE include)

install(DIRECTORY data/ DESTINATION data)

# This order is important as CMAKE needs to find certain files
add_subdirectory(lib)
add_subdirectory(src)

enable_testing()
add_subdirectory(tests)
