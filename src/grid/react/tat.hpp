#pragma once
#include "barrier.hpp"
#include "../../assertions.h"

namespace grid {
    namespace react {
        class tat_simple : public barrier{
            public:
                
            inline tat_simple(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s,std::shared_ptr<atom::Type> _to1, std::shared_ptr<atom::Type> _to2,double _maxdist,double _barrier,double _freq):barrier(f,s,_to1,_to2,_maxdist,_barrier,_freq) {};

            inline virtual double Chance(std::shared_ptr<grid::atom::Atom>& f, std::shared_ptr<grid::atom::Atom>& s, double distance) const {
                if(!AreAtomsOk(f,s)) return 0;
                
                double dE = getDEChecked(f,s);

                double E = -barrier_height+dE;

                double freq_elastic = freq * exp(-distance/Distance());

                double E_norm = E/(sgs::BOLZMAN*f->T());
                
                double freq_inelastic = freq*E_norm*1/(1-exp(-E_norm));

                double freqq = freq_elastic + freq_inelastic;

                assert_tst(freqq==freqq);

                return freqq;
            };
        };
    }
}