#pragma once

#include "../grid/grid.hpp"

namespace field {
    using grid::Grid;
    class Field {
        public:

        inline virtual void Apply(const Grid& g) const {};

        inline virtual void Apply(const geometry::Vector& pos,const std::shared_ptr<grid::atom::Atom>& patom) {};
    };
}