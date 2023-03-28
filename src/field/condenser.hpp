#pragma once

#include "field.hpp"

namespace field {
    using grid::Grid;
    class ZCondenser: public Field{
        double Uf;

        double start, end;
        public:
        inline ZCondenser(double _U, double _start, double _end): Uf(_U/(_end-_start)),start(_start),end(_end) {}

        inline void Apply(const Grid& g) const {
            auto field = this;
            g.for_each([field](const geometry::Vector& pos,std::shared_ptr<grid::atom::Atom> atom) {
                if(pos.z>=field->start && pos.z<=field->end) {
                    double dU = atom->U()+field->Uf*(pos.z-field->start);
                    atom->U(dU);
                }
            });
        };

        inline void Apply(const geometry::Vector& pos,const std::shared_ptr<grid::atom::Atom> patom) {
            if(pos.z>=start && pos.z<=end) 
                patom->U(patom->U()+Uf*(pos.z-start));

        };
    };
}