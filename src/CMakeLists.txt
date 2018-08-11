add_executable( stch
    ProgramManager.cpp
    main.cpp

    programs/ScorePose.cpp
    programs/AllScorePose.cpp
    programs/KBMinimize.cpp
    programs/PhysMinimize.cpp
)

target_include_directories(stch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(stch PUBLIC statchem pthread)
