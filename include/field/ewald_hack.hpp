#pragma once

#include <grid/grid.hpp>
#include <unordered_map>
#include <geometry/vector.hpp>

namespace field {
    using geometry::Vector;

    struct evald_chunk_data {
        public:
        double q;
    };
    class ewald_hack {
        std::unordered_map<Vector,evald_chunk_data> data;
        grid::Grid* grid;
        double chunk_size;
        bool cycle_x,cycle_y,cycle_z;


        public:

        inline ewald_hack(grid::Grid& g):grid(&g),data(),chunk_size(g.Chunk_Size()) {
            cycle_x = g.Cyclic<'x'>();
            cycle_y = g.Cyclic<'y'>();
            cycle_z = g.Cyclic<'z'>();
        };

        inline void calc_c() {
            for( auto&pchunk: grid->Chunks()) {
                auto& vchunk = pchunk.first;
                double q=0;
                for(auto& pa:*pchunk.second) {
                    auto& atom = pa.second;
                    q+=atom->Q();
                }


                data[vchunk] = {q};
            }
        };

        inline double calc_a(const Vector& pos) {
            double U = 0;
            auto cp = grid->calcChunkPos(pos);

            for(auto& [vec,cdata]:data) {
                
                if(vec!=cp){

                    auto p1 = vec+Vector(1,1,1)*chunk_size/2;
                    auto vc = grid->getMinDist(p1,pos);
                    double dist = vc.abs();
                    U += cdata.q/dist;


                }
                else {
                    auto& chunks = grid->Chunks();
                    auto pchunk = chunks.find(cp);
                    if(pchunk!=chunks.end()) {
                        for(auto& [vvec,atom]:*pchunk->second) {
                            if(vvec!=pos) {
                                auto vc = grid->getMinDist(vvec,pos);
                                U+=0.1*atom->Q()/(vc).abs();
                               
                            }
                        }
                    }
                }
            }
            return U;
        };

        inline void calc_all() {
            calc_c();
            for(auto iter = grid->begin();!(iter==grid->end());) {
                auto& vec = iter.aiter->first;
                auto& atom = iter.aiter->second;
                ++iter;
                atom->U(atom->U()+calc_a(vec));



            }        
        }
    };
}