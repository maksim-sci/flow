set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/geometry/vector.hpp"
                "${PROJECT_SOURCE_DIR}/include/sgs.hpp"
                "${PROJECT_SOURCE_DIR}/include/geometry/geometry.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/cubicchunk.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/atom/type.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/atom/atom.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/lattice.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/react/react.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/react/puasson.hpp"
                "${PROJECT_SOURCE_DIR}/include/grid/grid.hpp"
                "${PROJECT_SOURCE_DIR}/include/math/modulo.hpp"
                "${PROJECT_SOURCE_DIR}/include/field/condenser.hpp"
                "${PROJECT_SOURCE_DIR}/include/field/equal.hpp"
                "${PROJECT_SOURCE_DIR}/include/field/evald.hpp"
                "${PROJECT_SOURCE_DIR}/include/field/ewald_hack.hpp"
)

# Make an automatic library - will be static or dynamic based on user setting
add_library(flow_library
        geometry/vector.cpp
        geometry/geometry.cpp
        grid/cubicchunk.cpp
        grid/react.cpp
        grid/grid.cpp
        grid/lattice.cpp
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
