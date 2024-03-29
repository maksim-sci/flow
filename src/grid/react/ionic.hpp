#pragma once
#include "barrier.hpp"
#include "../../assertions.h"

namespace grid {
    namespace react {
        class ionic : public barrier{
            public:
                
            inline ionic(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s, std::shared_ptr<atom::Type> _to1, std::shared_ptr<atom::Type> _to2,double _maxdist,double _barrier,double _freq):barrier(f,s,_to1,_to2,_maxdist,_barrier,_freq) {};

            inline virtual double Chance(std::shared_ptr<grid::atom::Atom>& f, std::shared_ptr<grid::atom::Atom>& s, double distance) const {
                if(!AreAtomsOk(f,s)) return 0;
                
                double dE = getDEChecked(f,s);

                double E = -barrier_height+dE;
                double freqq = freq * exp( (E)/(sgs::BOLZMAN*f->T()) ) *exp(-distance/maxdist);

                assert_tst(freqq==freqq);

                return freqq;
            };
        };
    }
}