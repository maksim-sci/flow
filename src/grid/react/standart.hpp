#pragma once
#include "react.hpp"
#include "../../assertions.h"

namespace grid {
    namespace react {
        class standart : public react{
            public:
            double barrier;
            double freq;
            public:
                inline standart(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s,std::shared_ptr<atom::Type> _to1, std::shared_ptr<atom::Type> _to2,double _maxdist,double _barrier,double _freq):react(f,s,_to1,_to2,_maxdist),barrier(_barrier),freq(_freq) {};

                virtual double calcDE(const grid::atom::Atom* f, const grid::atom::Atom* s) const;

                inline virtual double getCutoff() const{
                    return sgs::ELVOLT*20;
                }

                virtual double getDEChecked(const grid::atom::Atom* f, const grid::atom::Atom* s) const;

                virtual double Chance(const grid::atom::Atom* f, const grid::atom::Atom* s, double distance) const;

                virtual void Apply(grid::atom::Atom* f, grid::atom::Atom* s) const ;
        };
    }
}
