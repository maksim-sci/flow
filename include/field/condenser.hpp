#pragma once

#include <grid/grid.hpp>

namespace field {
    using grid::Grid;
    class ZCondenser {
        double Uf;

        double start, end;
        public:
        inline ZCondenser(double _U, double _start, double _end): Uf(_U/(_end-_start)),start(_start),end(_end) {}

        inline void Apply(const Grid& g) const {
            auto iter = g.begin();
            auto iend = g.end();

            for (;!(iter==iend);) {
                auto& pos = iter.aiter->first;

                if(pos.z>=start && pos.z<=end) {
                    auto& patom = iter.aiter->second;
                    patom->U(patom->U()+Uf*(pos.z-start));
                }
                ++iter;
            }
        }
    };
}