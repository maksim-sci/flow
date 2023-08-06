#include "barrier.hpp"

namespace grid {
    namespace react {

        double barrier::calcDE(std::shared_ptr<grid::atom::Atom>& f,std::shared_ptr<grid::atom::Atom>& s) const{
            double dq = s->Q()-f->Q();
            double sq = s->Q();
            double fq = f->Q();
            double su = s->U();
            double fu = f->U();
            double dE = -(-fu*fq-su*sq+su*to2->Q()+fu*to1->Q())/2;

            return dE;
        }

        double barrier::getDEChecked(std::shared_ptr<grid::atom::Atom>& f,std::shared_ptr<grid::atom::Atom>& s) const {
            double dE = calcDE(f,s);

            return dE;
        }

        double barrier::Chance(std::shared_ptr<grid::atom::Atom>& f,std::shared_ptr<grid::atom::Atom>& s, double distance) const {
            if(!AreAtomsOk(f,s)) return 0;
            
            double dE = getDEChecked(f,s);

            double E = -barrier_height+dE;
            double freqq = freq * exp( (E)/(sgs::BOLZMAN*f->T()) );

            assert_tst(freqq==freqq);

            return freqq;
        };

        void barrier::Apply(std::shared_ptr<grid::atom::Atom>&f, std::shared_ptr<grid::atom::Atom>& s) const {
            if(!AreAtomsOk(f,s))return;
            f->Material(to1);
            s->Material(to2);
        };
    }
}