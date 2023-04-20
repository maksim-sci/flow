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

        public:

        ewald(double _real_cutoff, double _reciprocal_cutoff, int _calc_cnt, double _sigma,grid::Grid* _g):real_cutoff(_real_cutoff),reciprocal_cutoff(_reciprocal_cutoff),calc_cnt(_calc_cnt),sigma(_sigma),g(_g){};

        void apply();
        

    };
}