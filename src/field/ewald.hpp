#pragma once
#include "../grid/grid.hpp"
#include "../geometry/vector.hpp"

namespace field {

    class ewald {
        double real_cutoff;
        double reciprocal_cutoff;
        int calc_cnt;
        double kappa;
        double sigma;
        grid::Grid* g;

        double size_x;
        double size_y;
        double size_z;
        double recip_x;
        double recip_y;
        double recip_z;
        double V;
        geometry::Vector recip_size;

        public:

        ewald(double _real_cutoff, double _reciprocal_cutoff, int _calc_cnt, double _sigma,grid::Grid* _g):real_cutoff(_real_cutoff),reciprocal_cutoff(_reciprocal_cutoff),calc_cnt(_calc_cnt),sigma(_sigma),g(_g)
        {
        };

        void cache_grid_data() {
            const auto sizes = g->Sizes();
            const auto begin = g->Llim();
            size_x=sizes.x;
            size_y=sizes.y;
            size_z=sizes.z;
            V = sizes.x*sizes.y*sizes.z;
            recip_x = (2 * M_PI / V) * (size_y * size_z);
            recip_y = (2 * M_PI / V) * (size_z * size_x);
            recip_z = (2 * M_PI / V) * (size_x * size_y);

            recip_size=geometry::Vector{recip_x,recip_y,recip_z};
        };

        void apply() const;
        void add_charge(const geometry::Vector pos, double q) const;
        

    };
}