#pragma once
#include "../grid/grid.hpp"
#include "../geometry/vector.hpp"
#include <unordered_map>
namespace field {
    class ChunkProperties {
        double u;
    };
    class evald {
        double chunk_size;
        geometry::Vector translations;
        grid::Grid* grid;

        inline evald(grid::Grid& g):grid(&g),chunk_size(g.Chunk_Size()),translations({g.Chunk_Size(),g.Chunk_Size(),g.Chunk_Size()}){

        };

        inline void apply() {

        };

    };
}