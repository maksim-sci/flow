set(HEADER_LIST "geometry/vector.hpp"
                "sgs.hpp"
                "geometry/geometry.hpp"
                "grid/cubicchunk.hpp"
                "grid/atom/type.hpp"
                "grid/atom/atom.hpp"
                "grid/lattice.hpp"
                "grid/react/react.hpp"
                "grid/react/barrier.hpp"
                "grid/grid.hpp"
                "math/modulo.hpp"
                "field/condenser.hpp"
                "field/equal.hpp"
                "field/ewald_hack.hpp"
                "field/basic_electrode.hpp"
                "field/field.hpp"
                "field/summation.hpp"
                "field/ewald.hpp"
                "grid/react/ionic.hpp"
                "grid/react/tat.hpp"
                "grid/react/tat_simple.hpp"
                "algo/kmk.hpp"
)

# Make an automatic library - will be static or dynamic based on user setting
add_library(flow_library
        geometry/vector.cpp
        geometry/geometry.cpp
        grid/cubicchunk.cpp
        grid/react/react.cpp
        grid/react/barrier.cpp
        grid/grid.cpp
        grid/lattice.cpp
        field/ewald.cpp
        field/summation.cpp
        algo/kmk.cpp
        ${HEADER_LIST})


link_libraries(stdc++fs)

target_link_libraries(flow_library PRIVATE fmt::fmt)

# All users of this library will need at least C++17
target_compile_features(flow_library PUBLIC cxx_std_20)

target_include_directories(flow_library PUBLIC "./")