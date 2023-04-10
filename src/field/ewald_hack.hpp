#pragma once

#include "../grid/grid.hpp"
#include "../geometry/vector.hpp"
#include <unordered_map>

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

        inline void add_chunk_q(const geometry::Vector& pos,double dq) {
            if(dq==0) return;
            auto chunk = data.find(grid->calcChunkPos(pos));
            if(chunk!=data.end()) {
                chunk->second.q+=dq;
            }
        }

        inline void calc_chunk(const geometry::Vector& pos) {
            auto& chunks = grid->Chunks();
            auto pchunk = chunks.find(pos);
            if(pchunk!=chunks.end()) {
                auto& [vec,chunk] = *pchunk;
                for(auto iter = chunk->begin();iter!=chunk->end();iter++) {
                    auto& [atom_pos,atom] = *iter;
                    atom->U(calc_a(atom_pos));
                }
            }
        }

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

            for(auto& [vec,chunk_data]:data) {
             

                auto center = vec+Vector(1,1,1)*chunk_size/2;
                auto delta = grid->getMinDist(center,pos);
                double dist = delta.abs();
                U += chunk_data.q/dist;
            }
            return U;
        };

        inline void calc_all() {
            data.clear();
            data.reserve(10000);
            data.max_load_factor(0.25);
            calc_c();
            grid->for_each([this](const Vector& pos,std::shared_ptr<grid::atom::Atom> atom){
                double du = this->calc_a(pos);
                atom->U(atom->U()+du);
            });
        }
    };
}