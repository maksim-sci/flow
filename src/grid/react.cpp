
#include <grid/react/react.hpp>
#include <sgs.hpp>
namespace grid {
    namespace react {
        bool React::AreAtomsOk(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s)const {
            if(s==f) return false;
            if(f->Material()!=first) {return false;}
            if(s->Material()!=second) {return false;}

            return true;
        }

        double React::Chance(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s, double distance) const{
            if(!AreAtomsOk(f,s)) return 0;
            return 1*exp(s->U()-f->U()/sgs::BOLZMAN*s->T());
        }
    }
}