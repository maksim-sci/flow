#pragma once

#include "field.hpp"

namespace field {
    using grid::Grid;
    class basicElectrode: public Field{
        using geometry::Vector;

        double Uf;
        double EUf;
        Vector begin,ebd;
        public:
        inline basicElectrode(double _U, double _start, double _end,double _EUf,Vector  _begin, Vector _end): Uf(_U/(_end-_start)),start(_start),end(_end),EUf(_EUf/_end.z-_begin.z),begin(_begin),end(_end) {}

        inline void Apply(const Grid& g) const {
            auto iter = g.begin();
            auto iend = g.end();

            for (;!(iter==iend);) {
                auto& pos = iter.aiter->first;

                if(pos.z>=start && pos.z<=end) {
                    auto& patom = iter.aiter->second;
                    
                    double dU;
                    if(geometry::IsVectorInCube(patom,begin,end)) {
                        dU = patom->U()+EUf*(pos.z-begin.z);
                    }
                    else
                    {
                        dU = patom->U()+Uf*(pos.z-start);
                    }
                    patom->U(dU);
                    //fmt::print("result: {:e}\n",dU);
                }
                ++iter;
            }
        };

        inline void Apply(const geometry::Vector& pos,const std::shared_ptr<grid::atom::Atom>& patom) {
            if(pos.z>=start && pos.z<=end) {
                double dU;
                if(geometry::IsVectorInCube(pos,begin,end)) {
                    dU = patom->U()+EUf*(pos.z-begin.z);
                }
                else
                {
                    dU = patom->U()+Uf*(pos.z-start);
                }
                patom->U(dU);
            }
        };
    };
}