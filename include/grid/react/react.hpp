#pragma once
#include <memory>

#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
namespace grid {
    namespace react {

        class React {
        private:
            std::shared_ptr<atom::Type> first;
            std::shared_ptr<atom::Type> second;

            double maxdist;



        public:
            inline React(std::shared_ptr<atom::Type> f,std::shared_ptr<atom::Type> s,double _maxdist):first(f),second(s),maxdist(_maxdist) {};

            virtual void Apply(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s) const {};

            virtual bool AreAtomsOk(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s)const ;
            virtual double Chance(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s, double distance) const;
            inline double Distance() {
                return maxdist;
            };

        };
    }
}