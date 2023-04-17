
#include "react.hpp"
#include "../../sgs.hpp"
namespace grid {
    namespace react {
        bool react::AreAtomsOk(grid::atom::Atom* const f, grid::atom::Atom* const s)const {
            if(s==f) return false;
            if(f->Material()!=first) {return false;}
            if(s->Material()!=second) {return false;}

            return true;
        };

        double react::Chance(grid::atom::Atom* const f, grid::atom::Atom* const s, double distance) const{
            if(!AreAtomsOk(f,s)) return 0;
            return 1*exp(s->U()-f->U()/sgs::BOLZMAN*s->T());
        };
    }
}