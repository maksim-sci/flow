#pragma once
#include "react.hpp"

namespace grid {
    namespace react {
        class Puasson : public React{
            public:
            double barrier;
            double freq;
            public:
                inline Puasson(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s,std::shared_ptr<atom::Type> _to1, std::shared_ptr<atom::Type> _to2,double _maxdist,double _barrier,double _freq):React(f,s,_to1,_to2,_maxdist),barrier(_barrier),freq(_freq) {};

                inline virtual double Chance(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s, double distance) const {
                    if(!AreAtomsOk(f,s)) return 0;
                    if(distance>this->maxdist) return 0;
                    double dq = s->Q()-f->Q();
                    double sq = s->Q();
                    double fq = f->Q();
                    double su = s->U();
                    double fu = f->U();
                    double dE = fu*fq+su*sq-su*to2->Q()-fu*to1->Q();
                    dE*=-0.01;

                    if(dE>sgs::BOLZMAN*100) {
                        dE = sgs::BOLZMAN*100;
                    }

                    double E = -barrier+dE;
                    double freqq = freq * exp( (E)/(sgs::BOLZMAN*f->T()) );

                    if(freqq!=freqq) {
                        int* a = 0;
                        *a = 0;
                    }
                    return freqq;
                };

                inline virtual void Apply(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s) const {
                    if(!AreAtomsOk(f,s))return;
                    f->Material(to1);
                    s->Material(to2);
                };
        };
    }
}
