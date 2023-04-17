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
        g->for_each([&](const auto& pos, const auto& atom) mutable {
            data.push_back({atom,pos,0});
        });

        const auto sizes = g->Sizes();
        const auto begin = g->Llim();

        double V = sizes.x*sizes.y*sizes.z;

        Vector reciprocal_trans{1/sizes.x,1/sizes.y,1/sizes.z};
        reciprocal_trans*=2*M_PI;

        for(auto& atom_info_1:data) {
            auto& atom1 = (atom_info_1.atom);
            for(auto& atom_info_2:data) {
                auto& atom2 = (atom_info_2.atom);
                Vector delta = atom_info_2.real_pos-atom_info_1.real_pos;
                double u_reciprocal = 0;
                double u_real = 0;
                for(int x = -calc_cnt;x <= calc_cnt;x++) {
                    for(int y = -calc_cnt;y <= calc_cnt;y++) {
                        int z = 0;
                        //reciprocal space
                        {
                            Vector delta_rec {x*reciprocal_trans.x,y*reciprocal_trans.y,0};
                            double dist = delta_rec.abs();
                            if(x==0&&y==0){}
                            else {
                                if(dist<reciprocal_cutoff) {
                                    double diff = (atom1->Q()*atom2->Q())*(1/dist)*exp(-dist/(4.0*kappa*kappa))*cos(delta_rec.scalar_mul(delta));
                                    assert_tst(diff==diff);
                                    u_reciprocal+= diff;
                                }
                            }
                        }
                        

                        //real space
                        {
                            Vector delta_real {x*sizes.x,y*sizes.y,0};
                            double dist = delta_real.abs();

                            double len = (delta+delta_real).abs();
                            if(dist<real_cutoff) {
                                if(atom1!=atom2) {
                                    double diff = (atom1->Q()*atom2->Q()) * erfc(kappa*len)/len;
                                    assert_tst(diff==diff);
                                    u_real += diff;
                                }
                            }
                        }
                        
                    }
                }
                double u = u_reciprocal*4.0*M_PI * 0.5/V+u_real/2;

                atom2->U(u/sgs::ELCHARGE);
            }
        }

    };
}
