if (WIN32)
    if (NOT ${STCH_BUILD_STATIC_EXECUTABLE})
            add_definitions( /DSTCH_SHARED_LIBRARY)
    else()
            add_definitions( /DSTCH_STATIC_LIBRARY)
    endif()
endif(WIN32)

file(GLOB_RECURSE STCH_LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/**.cpp)

if (NOT ${STCH_BUILD_STATIC_EXECUTABLE})

    add_library(
        statchem SHARED
            ${STCH_LIB_SOURCES}
    )

    target_link_libraries ( statchem
        ${OPENMM_LIBRARY}
        ${Boost_LIBRARIES}
        ${GSL_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
    )

else()
    add_library(
        statchem STATIC
            ${STCH_LIB_SOURCES}
    )
endif()

install(
    TARGETS
        statchem
    LIBRARY DESTINATION
        ${CMAKE_INSTALL_PREFIX}/lib
    ARCHIVE DESTINATION
        ${CMAKE_INSTALL_PREFIX}/lib
    RUNTIME DESTINATION
        ${CMAKE_INSTALL_PREFIX}/bin
)

if (${USING_INTERNAL_GSL})
    add_dependencies(statchem GSL)
endif()

if (${USING_INTERNAL_OPENMM})
    add_dependencies(statchem OpenMM_Build)
endif()
