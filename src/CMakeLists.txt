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
)

target_include_directories(stch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(stch PUBLIC statchem pthread)
