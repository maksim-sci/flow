#pragma once
#include <memory>
#include <string>

#include "../atom/type.hpp"
#include "../atom/atom.hpp"
namespace grid {
    namespace react {

        class react {
        private:
            public:
            std::shared_ptr<atom::Type> first;
            std::shared_ptr<atom::Type> second;

            
            std::shared_ptr<grid::atom::Type> to1;
            std::shared_ptr<grid::atom::Type> to2;

            std::string name;

            double maxdist;



        public:
            inline react(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s,std::shared_ptr<atom::Type> _to1, std::shared_ptr<atom::Type> _to2,double _maxdist):name(""),first(f),second(s),maxdist(_maxdist),to1(_to1),to2(_to2) {};

            virtual void Apply(const grid::atom::Atom* f, const grid::atom::Atom* s) const {};

            virtual bool AreAtomsOk(grid::atom::Atom* const f, grid::atom::Atom* const s)const ;
            virtual double Chance(grid::atom::Atom* const f, grid::atom::Atom* const s, double distance) const;
            inline double Distance() {
                return maxdist;
            };
            virtual inline void Name(std::string& n) {name = n;};
            virtual inline const std::string& Name() {return name;};

        };
    }
}