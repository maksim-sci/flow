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

            inline virtual void Apply(grid::atom::Atom* f, grid::atom::Atom* s) const {};
            //inline virtual void Apply(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s) const {Apply(f.get(),s.get());};


            virtual bool AreAtomsOk(const grid::atom::Atom* f, const grid::atom::Atom* s)const ;
            //inline virtual bool AreAtomsOk(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s)const {return AreAtomsOk(f.get(),s.get());};
            virtual double Chance(const grid::atom::Atom* f, const grid::atom::Atom* s, double distance) const;
            //inline virtual double Chance(std::shared_ptr<grid::atom::Atom> f, std::shared_ptr<grid::atom::Atom> s, double distance) const {return Chance(f.get(),s.get(),distance);};
            inline double Distance() {
                return maxdist;
            };
            virtual inline void Name(std::string& n) {name = n;};
            virtual inline const std::string& Name() {return name;};

        };
    }
}