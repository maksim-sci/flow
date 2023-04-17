#include "standart.hpp"

namespace grid {
    namespace react {

        double standart::calcDE(grid::atom::Atom* const f,  grid::atom::Atom* const s) const{
            double dq = s->Q()-f->Q();
            double sq = s->Q();
            double fq = f->Q();
            double su = s->U();
            double fu = f->U();
            double dE = (-fu*fq-su*sq+su*to2->Q()+fu*to1->Q()*sgs::ELVOLT/(sgs::VOLT*sgs::ELCHARGE))/2;

            return dE;
        }

        double standart::getDEChecked(grid::atom::Atom* const f, grid::atom::Atom* const s) const {
            double dE = calcDE(f,s);

            return dE;
        }

        double standart::Chance(grid::atom::Atom* const f, grid::atom::Atom* const s, double distance) const {
            if(!AreAtomsOk(f,s)) return 0;
            
            double dE = getDEChecked(f,s);

            double E = -barrier+dE;
            double freqq = freq * exp( (E)/(sgs::BOLZMAN*f->T()) );

            assert_tst(freqq==freqq);

            return freqq;
        };

        void standart::Apply(grid::atom::Atom* f, grid::atom::Atom* s) const {
            if(!AreAtomsOk(f,s))return;
            f->Material(to1);
            s->Material(to2);
        };
    }
}