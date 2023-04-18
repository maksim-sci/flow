#pragma once

#include "../grid/grid.hpp"

namespace field {
    using grid::Grid;
    class Equal {
        double U;

        double start, end;
        public:
        inline Equal(double _U, double _start, double _end):U(_U),start(_start),end(_end) {};
        inline void Apply(const Grid& g) const
        {
            g.for_each([this](const geometry::Vector& pos,std::shared_ptr<grid::atom::Atom>& atom)
            {
                if(pos.z>=this->start && pos.z<=this->end) {
                    atom->U(this->U);
                }
            });
        };

        inline void Apply(const geometry::Vector& pos,const std::shared_ptr<grid::atom::Atom>& patom) {
            if(pos.z>=start && pos.z<=end) 
                patom->U(U);
        };
    };
}