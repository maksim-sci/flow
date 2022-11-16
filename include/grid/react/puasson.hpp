#pragma once
#include "react.hpp"

namespace grid {
    namespace react {
        class Puasson : public React{
            double barrier;
            double freq;
            std::shared_ptr<grid::atom::Type> to1;
            std::shared_ptr<grid::atom::Type> to2;
            public:
                inline Puasson(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s,std::shared_ptr<atom::Type> _to1, std::shared_ptr<atom::Type> _to2,double _maxdist,double _barrier,double _freq):React(f,s,_maxdist),barrier(_barrier),to1(_to1),to2(_to2),freq(_freq) {};
                inline virtual double Chance(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s, double distance) const {
                    if(!AreAtomsOk(f,s)) return 0;
                    double E = barrier+s->U()-f->U();
                    return freq*(E)/((sgs::PLANCK)*(1-exp(-E/(sgs::BOLZMAN*f->T()))));
                };
                inline virtual void Apply(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s) const {
                    f->Material(to1);
                    s->Material(to2);
                };
        };
    }
}
