


#include <iostream>
#include <sgs.hpp>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <functional>
#include <memory>

#include <fmt/format.h>

#include <geometry/vector.hpp>
#include <geometry/geometry.hpp>
#include <grid/cubicchunk.hpp>
#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/lattice.hpp>
#include <grid/react/react.hpp>
#include <grid/grid.hpp>
#include <sgs.hpp>
#include <math/modulo.hpp>
#include <field/ewald_hack.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>
#include <filesystem>


namespace fs = std::filesystem;

using std::string;
using geometry::Vector;
using geometry::Geometry;
using grid::chunk::CubicChunk;
using grid::atom::Atom;
using grid::atom::Type;
using grid::Lattice;
using grid::atom::Atom;
using grid::react::React;
using grid::Grid;

#define assert_eq(a,b) {if(!((a)==(b))) {fmt::print("{}:{} {}!={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));exit(-1);}}
#define assert_neq(a,b) {if((a)==(b)) {fmt::print("{}:{} {}=={}: {}!={}\n",__FILE__,__LINE__,#a,#b,(a),(b));exit(-1);}}
#define assert_tst(a) {if(!(a)) {fmt::print("{}:{} (!a)=true: {}",__FILE__,__LINE__,#a); exit(-1);}}
void printvoltage(fs::path pout,Grid& g)
{
        std::ofstream out(pout);
        out << "x y z c\n";
        for (auto iter = g.begin(); !(iter == g.end());)
        {
            auto &[vec, a] = *iter.aiter;
            auto &[x, y, z] = vec;
            out << fmt::format("{} {} {} {}\n", x, y, z, a->U());
            ++iter;
        }
        out.close();
};

double calc_Energy(Grid& g) {
    double E = 0;
    for (auto iter = g.begin(); !(iter == g.end());)
    {
        auto &[vec, a] = *iter.aiter;
        E+=a->U()*a->Q();
        
        ++iter;
    }

    return E;
}


int main() {
    auto TElectrode = std::make_shared<Type>(0, __COUNTER__, "El");
    auto Oxygen = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "O");
    auto Oxygen_Intersittal = std::make_shared<Type>(-2 * sgs::ELCHARGE, __COUNTER__, "OI");
    auto Hafnium = std::make_shared<Type>(+2 * sgs::ELCHARGE, __COUNTER__, "Hf");
    auto OxygenVacancy_Neutral = std::make_shared<Type>(0, __COUNTER__, "Vo");
    auto OxygenVacancy_Neutral_Intersitial = std::make_shared<Type>(0, __COUNTER__, "Vi");
    auto OxygenVacancy_Charged = std::make_shared<Type>(+2 * sgs::ELCHARGE, __COUNTER__, "Vp");
    auto ElectrodePositive = std::make_shared<Type>(+1 * sgs::ELCHARGE, __COUNTER__, "Ep");
    auto ElectrodeNegative = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "Em");



    Vector HfO2A(5.069186, 0.000000, -0.864173);
    Vector HfO2B(0.000000, 5.195148, 0.000000);
    Vector HfO2C(0.000000, 0.000000, 5.326038);

    HfO2A *= sgs::ANGSTROM;
    HfO2B *= sgs::ANGSTROM;
    HfO2C *= sgs::ANGSTROM;

    Geometry geoHfO2(HfO2A, HfO2B, HfO2C);
    Lattice lHfO2(geoHfO2);

    lHfO2.add({0.724041, 0.457319, 0.292109}, Hafnium);
    lHfO2.add({0.275959, 0.957319, 0.207891}, Hafnium);
    lHfO2.add({0.275959, 0.542681, 0.707891}, Hafnium);
    lHfO2.add({0.724041, 0.042681, 0.792109}, Hafnium);
    lHfO2.add({0.551113, 0.742603, 0.022292}, Oxygen);
    lHfO2.add({0.448887, 0.242603, 0.477708}, Oxygen);
    lHfO2.add({0.448887, 0.257397, 0.977708}, Oxygen);
    lHfO2.add({0.551113, 0.757397, 0.522292}, Oxygen);
    lHfO2.add({0.932151, 0.330122, 0.652906}, Oxygen);
    lHfO2.add({0.067849, 0.830122, 0.847094}, Oxygen);
    lHfO2.add({0.067849, 0.669878, 0.347094}, Oxygen);
    lHfO2.add({0.932151, 0.169878, 0.152906}, Oxygen);


    Grid g(sgs::ANGSTROM*10);

    g.AddLattice({0,0,0},sgs::ANGSTROM*(Vector{5.05,5.18,5.3}),lHfO2);

    g.Cyclic<'x'> (true);
    g.Cyclic<'y'> (true);
    g.Cyclic<'z'> (true);




    field::ewald_hack hh(g);
    hh.calc_all();
    std::cout<<"atoms: "<<g.count()<<std::endl;
    std::cout<<"calculated energy: "<<calc_Energy(g)/sgs::ELVOLT<<"  Ev"<<std::endl;

    printvoltage("volt.txt",g);


}