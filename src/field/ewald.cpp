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
        std::vector<evald_atom_info> data;

        const auto sizes = g->Sizes();
        const auto begin = g->Llim();
        const auto& [size_x,size_y,size_z] = sizes;
        double V = sizes.x*sizes.y*sizes.z;

        constexpr double sqrt_2 = std::sqrt(2);
        constexpr double min_dist = sgs::ANGSTROM*0.1;

        g->for_each([&](const auto& pos1, const auto& atom1) mutable {
            double u_nn = 0;

            //real space
            g->for_each(pos1,real_cutoff,[&](const auto& pos2, const auto& atom2) mutable {
                if(atom1==atom2) return;
                double delta = (pos1-pos2).abs();
                if(delta<real_cutoff && delta > min_dist) {
                    u_nn+=atom2->Q()/delta*erf(delta/(sqrt_2*sigma));
                }
            });

            //reciprocal space
            g->for_each([&](const auto& pos2, const auto& atom2) mutable {
                for(int x = -calc_cnt;x <= calc_cnt;x++) {
                    for(int y = -calc_cnt;y <= calc_cnt;y++) {
                        if(x==0&&y==0&&atom1==atom2) return;
                        constexpr int z = 0;
                        double delta = (pos1-pos2+Vector{x*size_x,y*size_y,z*size_z}).abs();
                        if(delta<reciprocal_cutoff && delta > min_dist) {
                            u_nn+=atom2->Q()/delta*erfc(delta/(sqrt_2*sigma));
                        }
                    }
                }
            });

            atom1->U(u_nn/(4*M_PI));
        });

    };
}
