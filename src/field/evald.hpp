#pragma once
#include "../grid/grid.hpp"
#include "../geometry/vector.hpp"

namespace field {

    class evald {
        double real_cutoff;
        double reciprocal_cutoff;
        int calc_cnt;
        double kappa;
        grid::Grid* g;


        evald(double _real_cutoff, double _reciprocal_cutoff, int _calc_cnt, double _kappa,grid::Grid* _g):real_cutoff(_real_cutoff),reciprocal_cutoff(_reciprocal_cutoff),calc_cnt(_calc_cnt),kappa(_kappa),g(_g){};

        void apply();
        

    };
}