add_executable( stch
    ProgramManager.cpp
    main.cpp

    programs/ScorePose.cpp
    programs/AllScorePose.cpp
    programs/KBMinimize.cpp
    programs/KBDynamics.cpp
    programs/PhysMinimize.cpp
    programs/PhysDynamics.cpp
    programs/MakeObjective.cpp
    programs/AssignAtomTypes.cpp
)

target_include_directories(stch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (MSVC)
    target_link_libraries(stch PUBLIC statchem)
else()
    target_link_libraries(stch PUBLIC statchem pthread)
endif()
install(TARGETS stch DESTINATION bin)
