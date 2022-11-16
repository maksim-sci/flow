#pragma once

#include <grid/grid.hpp>

namespace field {
    using grid::Grid;
    class Equal {
        double U;

        double start, end;
        public:
        inline Equal(double _U, double _start, double _end):U(_U),start(_start),end(_end) {};
        inline void Apply(const Grid& g) const {
            auto iter = g.begin();
            auto iend = g.end();

            for (;!(iter==iend);) {
                auto& pos = iter.aiter->first;

                if(pos.z>=start && pos.z<=end) {
                    auto& patom = iter.aiter->second;
                    patom->U(U);
                }
                ++iter;
            }
        }
    };
}