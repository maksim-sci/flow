set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/section.hpp"
                "${PROJECT_SOURCE_DIR}/include/atom.hpp"
                "${PROJECT_SOURCE_DIR}/include/params.hpp"
                "${PROJECT_SOURCE_DIR}/include/point.hpp"
                "${PROJECT_SOURCE_DIR}/include/types/atom.hpp"
                "${PROJECT_SOURCE_DIR}/include/types/reaction.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid.hpp"
                "${PROJECT_SOURCE_DIR}/include/flow.hpp"
                "${PROJECT_SOURCE_DIR}/include/tension.hpp"
                "${PROJECT_SOURCE_DIR}/include/app.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(flow_library
        section.cpp
        type_atom.cpp
        type_reaction.cpp
        point.cpp
        atom.cpp
        params.cpp
        grid.cpp
        flow.cpp
        tension.cpp
        app.cpp
        ${HEADER_LIST})


link_libraries(stdc++fs)

target_link_libraries(flow_library PRIVATE fmt::fmt)

# We need this directory, and users of our library will need it too
target_include_directories(flow_library PUBLIC ../include)

# All users of this library will need at least C++17
target_compile_features(flow_library PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
source_group(
        TREE "${PROJECT_SOURCE_DIR}/include"
        PREFIX "Header Files"
        FILES ${HEADER_LIST})
