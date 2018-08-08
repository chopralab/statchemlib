include(ExternalProject)

set(USING_INTERNAL_GSL on CACHE BOOL "Using Internal GSL")
mark_as_advanced(USING_INTERNAL_GSL)

ExternalProject_Add( GSL
    URL https://github.com/ampl/gsl/archive/644e768630841bd085cb7121085a688c4ff424d0.tar.gz
    URL_HASH SHA1=de333ea777c1880edb9a59d68c8ed3e8cc91810c
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsl
    CMAKE_CACHE_ARGS    -DGSL_DISABLE_TESTS:Bool=on
                        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/usr/
                        -DGSL_DISABLE_WARNINGS:Bool=on
                        -DBUILD_SHARED_LIBS:Bool=on
                        -DCMAKE_BUILD_TYPE:String=Release
)

add_library(gsl SHARED IMPORTED)
add_library(gslcblas SHARED IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/usr/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gsl${CMAKE_SHARED_LIBRARY_SUFFIX})
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/usr/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gslcblas${CMAKE_SHARED_LIBRARY_SUFFIX})
add_dependencies(gsl GSL)
add_dependencies(gslcblas GSL)

set(GSL_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/usr/include)

if (WIN32)

    set(GSL_LIBRARIES ${CMAKE_BINARY_DIR}/usr/bin/${CMAKE_SHARED_LIBRARY_PREFIX}gsl${CMAKE_SHARED_LIBRARY_SUFFIX})

else()

    set(GSL_LIBRARIES ${CMAKE_BINARY_DIR}/usr/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gsl${CMAKE_SHARED_LIBRARY_SUFFIX}
                      ${CMAKE_BINARY_DIR}/usr/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gslcblas${CMAKE_SHARED_LIBRARY_SUFFIX}
)

endif()

mark_as_advanced(GSL_INCLUDE_DIRS)
mark_as_advanced(GSL_LIBRARIES)

install(
    FILES
        ${GSL_LIBRARIES}
    DESTINATION
        ${CMAKE_INSTALL_PREFIX}/lib
)
