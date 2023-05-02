#include "ewald.hpp"
#include "../sgs.hpp"

#include <unordered_map>
#include "../m_pi.h"

#include "../assertions.h"


using geometry::Vector;

namespace field {

    struct evald_atom_info {
        std::shared_ptr<grid::atom::Atom> atom;
        geometry::Vector real_pos;
        double u;
    };

    void ewald::apply() {

        const auto sizes = g->Sizes();
        const auto begin = g->Llim();
        const auto& [size_x,size_y,size_z] = sizes;
        double V = sizes.x*sizes.y*sizes.z;

        double recip_x,recip_y,recip_z;
        recip_x = (2 * M_PI / V) * (size_y * size_z);
        recip_y = (2 * M_PI / V) * (size_z * size_x);
        recip_z = (2 * M_PI / V) * (size_x * size_y);

        Vector recip_size{recip_x,recip_y,recip_z};

        constexpr double cutoff = sgs::ANGSTROM*0.5;

        g->for_each([&](const auto& pos1, const auto& atom1) mutable {
            double u1 = 0;
            double u2 = 0;


            //real space
            g->for_each(pos1,real_cutoff,[&](const auto& pos2, const auto& atom2) mutable {
                Vector delta = pos2-pos1;
                double dist = delta.abs();
                if(dist>cutoff) {
                    constexpr double pi_4 = 1/(4*M_PI);
                    constexpr double sqrt_2 = std::sqrt(2);
                    u1+=pi_4*atom2->Q()/dist*erfc(dist/(sqrt_2*sigma));
                }
            });

            //reciprocal space
            g->for_each([&](const auto& pos2, const auto& atom2) mutable {
                Vector delta = pos2-pos1;

                for(int i = -calc_cnt;i<calc_cnt;i++) {
                    for(int j = -calc_cnt;j<calc_cnt;j++) {
                        constexpr int z = 0;

                        if(i==0&j==0) {}
                        else{
                            Vector delta_translated = pos2 + Vector{i*size_x,j*size_y,0};

                            Vector reciprocal_delta = Vector{recip_x*size_x/delta_translated.x,recip_y*size_y/delta_translated.y,recip_z*size_z/delta_translated.z};

                            double reciprocal_delta_abs = reciprocal_delta.abs();

                            if(reciprocal_delta_abs>1e+20) continue;


                            double reciprocal_delta_sqr = reciprocal_delta_abs*reciprocal_delta_abs;
                            double phase = reciprocal_delta.scalar_mul(delta);
                            double exp_factor = -sigma*sigma*reciprocal_delta_sqr/2;
                            u2+=(
                                atom2->Q()/reciprocal_delta_sqr*
                                cos(phase)*
                                exp(exp_factor)
                            );
                        }
                    }
                }
            });

            u2/=V;

            double u = u1+u2;

            atom1->U(u);
        });

    };
}
