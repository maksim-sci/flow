set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/section/section.hpp"
                "${PROJECT_SOURCE_DIR}/include/atom/atom.hpp"
                "${PROJECT_SOURCE_DIR}/include/params/params.hpp"
                "${PROJECT_SOURCE_DIR}/include/point/point.hpp"
                "${PROJECT_SOURCE_DIR}/include/types/atom.hpp"
                "${PROJECT_SOURCE_DIR}/include/types/reaction.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/grid.hpp"
                "${PROJECT_SOURCE_DIR}/include/flow/flow.hpp"
                "${PROJECT_SOURCE_DIR}/include/tension/tension.hpp"
                "${PROJECT_SOURCE_DIR}/include/app/app.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(vovka_library
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

target_link_libraries(vovka_library PRIVATE fmt::fmt)

# We need this directory, and users of our library will need it too
target_include_directories(vovka_library PUBLIC ../include)

# All users of this library will need at least C++17
target_compile_features(vovka_library PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
source_group(
        TREE "${PROJECT_SOURCE_DIR}/include"
        PREFIX "Header Files"
        FILES ${HEADER_LIST})
