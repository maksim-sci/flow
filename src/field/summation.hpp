
#pragma once
#include "../grid/grid.hpp"
#include "../geometry/vector.hpp"

namespace field {

    class summation {
        grid::Grid* g;

        public:

        summation(grid::Grid* _g):g(_g){};

        void apply();
        void add_charge(geometry::Vector pos, double q);

    };
}