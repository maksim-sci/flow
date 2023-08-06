#include "ewald.hpp"
#include "../sgs.hpp"

#include <unordered_map>
#include "../m_pi.h"

#include "../assertions.h"


using geometry::Vector;


namespace field {

    constexpr double cutoff = sgs::ANGSTROM*0.5;

    void ewald::add_charge(const geometry::Vector pos, double q) const{

        g->for_each([&](const auto& pos1, const auto& atom) mutable {
            double u1 = 0;
            double u2 = 0;


            Vector delta = pos-pos1;

            //real space
            double dist = delta.abs();
            if(dist>cutoff) {
                constexpr double pi_4 = 1/(4*M_PI);
                constexpr double sqrt_2 = std::sqrt(2);
                u1+=pi_4*q/dist*erfc(dist/(sqrt_2*sigma));
            }

            //reciprocal space

            for(int i = -calc_cnt;i<calc_cnt;i++) {
                for(int j = -calc_cnt;j<calc_cnt;j++) {
                    constexpr int z = 0;

                    if(i==0&j==0) {}
                    else{
                        Vector delta_translated = pos + Vector{i*size_x,j*size_y,0};

                        Vector reciprocal_delta = Vector{recip_x*size_x/delta_translated.x,recip_y*size_y/delta_translated.y,recip_z*size_z/delta_translated.z};

                        double reciprocal_delta_abs = reciprocal_delta.abs();

                        if(reciprocal_delta_abs>1e+20) continue;


                        double reciprocal_delta_sqr = reciprocal_delta_abs*reciprocal_delta_abs;
                        double phase = reciprocal_delta.scalar_mul(delta);
                        double exp_factor = -sigma*sigma*reciprocal_delta_sqr/2;
                        u2+=(
                            q/reciprocal_delta_sqr*
                            cos(phase)*
                            exp(exp_factor)
                        );
                    }
                }
            }


            atom->U(atom->U()+u1+u2);
        });


    }

    void ewald::apply() const{

        g->for_each([&](const auto& pos1, const auto& atom) mutable {
            add_charge(pos1,atom->Q());
        });

    };
}
