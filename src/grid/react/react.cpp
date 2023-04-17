
#include "react.hpp"
#include "../../sgs.hpp"
namespace grid {
    namespace react {
        bool react::AreAtomsOk(const grid::atom::Atom* f,const  grid::atom::Atom* s)const {
            if(s==f) return false;
            if(f->Material()!=first) {return false;}
            if(s->Material()!=second) {return false;}

            return true;
        };

        double react::Chance(const grid::atom::Atom* f,const  grid::atom::Atom* s, double distance) const{
            if(!AreAtomsOk(f,s)) return 0;
            return 1*exp(s->U()-f->U()/sgs::BOLZMAN*s->T());
        };
    }
}