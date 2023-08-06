#include "ewald.hpp"
#include "../sgs.hpp"

#include <unordered_map>
#include "../m_pi.h"

#include "../assertions.h"

#include "summation.hpp"


using geometry::Vector;


namespace field {
    void summation::apply() {
        g->for_each([&](const auto& pos, const auto& atom)mutable {
            add_charge(pos,atom->Q());
        });

    };
    void summation::add_charge(geometry::Vector pos, double q) {
        g->for_each([&](const auto& pos1,const auto& atom) mutable {
            constexpr double eps = 5e-2;
            auto [Lx,Ly,Lz] = g->Sizes();

            auto calc_row = [&](double y)->double {
                Vector post = Vector{0,y*Ly,0}+pos;
                Vector a = pos1-post;
                double t = a.y*a.y+a.z*a.z;


                double sum = 0;
                double x = 1;

                double d1,d2;
                do {

                    auto pair_calc = [&](double x)->double {
                        double ap = a.x-Lx*x;
                        double am = a.x+Lx*x;

                        double sp = sqrt(ap*ap+t);
                        double sm = sqrt(am*am+t);
                        
                        return (sp+sm)/(sp*sm);
                    };

                    d1 = pair_calc(x);
                    d2 = pair_calc(x+1);

                    x+=2;

                    sum+=d1;
                    sum+=d2;


                } while((d2-d1)/d2>eps);

                
                sum=q*sum+((a.abs()>sgs::ANGSTROM*0.01)?q/a.abs():0);

                return sum;
            };

            auto calc_row_pair = [&](double y)->double {
                return calc_row(y)+calc_row(-y);
            };

            double y = 1;

            double sum = calc_row(0);

            double d1,d2;
            do {
                d1 = calc_row_pair(y);
                d2 = calc_row_pair(y+1);
                sum+=d1;
                sum+=d2;

                y+=2;

            } while((d1-d2)/d2>eps);

            atom->U(atom->U()+sum);
        });
    };
}
