

#include <iostream>


#include <sgs.hpp>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <functional>
#include <memory>
#include <filesystem>

#include <fmt/format.h>

#include <geometry/vector.hpp>
#include <geometry/geometry.hpp>
#include <grid/cubicchunk.hpp>
#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/lattice.hpp>
#include <grid/react/react.hpp>
#include <grid/react/puasson.hpp>
#include <grid/grid.hpp>
#include <sgs.hpp>
#include <math/modulo.hpp>
#include <field/ewald_hack.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>

using geometry::Vector;
using geometry::Geometry;
using grid::atom::Atom;
using grid::atom::Type;
using std::string;
using grid::Grid;
using grid::Lattice;

struct kmk_data {
    std::shared_ptr<grid::react::React> react;
    std::shared_ptr<grid::atom::Atom> f;
    std::shared_ptr<grid::atom::Atom> s;
    double chance;
};

void grid_like_final_ex() {
    Grid g(sgs::ANGSTROM*10);


    auto TElectrode = std::make_shared<Type>(0,__COUNTER__,"El");
    auto Oxygen = std::make_shared<Type>(0,__COUNTER__,"O");
    auto Oxygen_Intersittal = std::make_shared<Type>(-2*sgs::ELCHARGE,__COUNTER__,"O");
    auto Hafnium = std::make_shared<Type>(0,__COUNTER__,"Hf");
    auto OxygenVacancy_Neutral = std::make_shared<Type>(0,__COUNTER__,"Vo");
    auto OxygenVacancy_Charged = std::make_shared<Type>(+2*sgs::ELCHARGE,__COUNTER__,"Vo");
    auto ElectrodePositive = std::make_shared<Type>(+1*sgs::ELCHARGE,__COUNTER__,"Ep");
    auto ElectrodeNegative = std::make_shared<Type>(-1*sgs::ELCHARGE,__COUNTER__,"Em");

    double ela = 2.866*sgs::ANGSTROM;

    Geometry geoElectrode(Vector(ela,0,0),Vector(0,ela,0),Vector(0,0,ela));

    Lattice lElectrode(geoElectrode);
    lElectrode.add(Vector(0,0,0),TElectrode);

    Vector HfO2A(5.069186, 0.000000, -0.864173);
    Vector HfO2B(0.000000, 5.195148, 0.000000);
    Vector HfO2C(0.000000, 0.000000, 5.326038);

    HfO2A*=sgs::ANGSTROM;
    HfO2B*=sgs::ANGSTROM;
    HfO2C*=sgs::ANGSTROM;

    Geometry geoHfO2(HfO2A,HfO2B,HfO2C);
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
    lHfO2.add({0.522616, 0.110041, 0.782365}, OxygenVacancy_Neutral);
    lHfO2.add({0.500000, 0.639959, 0.250000}, OxygenVacancy_Neutral);
    lHfO2.add({0.833333, 0.000000, 0.500000}, OxygenVacancy_Neutral);
    lHfO2.add({0.747211, 0.514587, 0.862533}, OxygenVacancy_Neutral);
    lHfO2.add({0.666667, 0.000000, 0.166667}, OxygenVacancy_Neutral);
    lHfO2.add({0.086123, 0.152080, 0.470800}, OxygenVacancy_Neutral);
    lHfO2.add({0.283925, 0.625000, 0.048547}, OxygenVacancy_Neutral);
    lHfO2.add({0.904740, 0.646880, 0.592346}, OxygenVacancy_Neutral);


    double size_x = 30*sgs::ANGSTROM;
    double size_y = 30*sgs::ANGSTROM;

    double electrode_z = 5*sgs::ANGSTROM;

    double oxyde_z = 20*sgs::ANGSTROM+electrode_z;

    double electrode2_z = 5*sgs::ANGSTROM+oxyde_z;


    g.AddLattice({0.,0,0},Vector(size_x,size_y,electrode_z),lElectrode);
    g.AddLattice({0,0,electrode_z},Vector(size_x,size_y,oxyde_z),lHfO2);
    g.AddLattice({0,0,oxyde_z},Vector(size_x,size_y,electrode2_z),lElectrode);

    double size_electrode = sgs::ANGSTROM*10;
    double size_electrode_z = sgs::ANGSTROM*5;
    double electrode_begin_xy = sgs::ANGSTROM*7;
    double electrode_begin_z = electrode_z;

    Vector ElCube_s(electrode_begin_xy,electrode_begin_xy,electrode_begin_z);
    Vector ElCube_e(electrode_begin_xy+size_electrode,electrode_begin_xy+size_electrode,electrode_begin_z+size_electrode_z);

    g.ClearParallelep(ElCube_s,ElCube_e);

    g.AddLattice(ElCube_s,ElCube_e,lElectrode);

    grid::react::Puasson R1(Oxygen,OxygenVacancy_Neutral,OxygenVacancy_Charged,Oxygen_Intersittal,5*sgs::ANGSTROM,sgs::ELVOLT*0.83,1e+13);
    grid::react::Puasson R2(OxygenVacancy_Neutral,Oxygen_Intersittal,Oxygen_Intersittal,OxygenVacancy_Neutral,5*sgs::ANGSTROM,sgs::ELVOLT*0.8,1e+13);
    grid::react::Puasson R3(OxygenVacancy_Charged,Oxygen_Intersittal,Oxygen,OxygenVacancy_Neutral,5*sgs::ANGSTROM,sgs::ELVOLT*0.83,1e+13);
    grid::react::Puasson R4(Oxygen_Intersittal,OxygenVacancy_Neutral,OxygenVacancy_Neutral,Oxygen_Intersittal,5*sgs::ANGSTROM,sgs::ELVOLT*0.8,1e+13);

    
    grid::react::Puasson E1(OxygenVacancy_Charged,OxygenVacancy_Neutral,OxygenVacancy_Neutral,OxygenVacancy_Charged,5*sgs::ANGSTROM,sgs::ELVOLT*0.5,1e+15);
    grid::react::Puasson E2(OxygenVacancy_Neutral,ElectrodeNegative,OxygenVacancy_Charged,ElectrodeNegative,5*sgs::ANGSTROM,sgs::ELVOLT*0.5,1e+15);
    grid::react::Puasson E3(OxygenVacancy_Charged,ElectrodePositive,OxygenVacancy_Neutral,ElectrodePositive,5*sgs::ANGSTROM,sgs::ELVOLT*0.5,1e+15);

    namespace fs = std::filesystem;
    
    auto outfolder = fs::absolute("C:/Users/sci/Desktop/modeltest/test");


    fs::create_directories(outfolder);

    fs::path out(outfolder);
    auto firstout = out;
    firstout/="first.xyz";

    auto printgrid = [](const Grid& g,fs::path pout){
        std::ofstream out(pout);

        out<<g.to_xyz(1/sgs::ANGSTROM);

        out.close();
    };

    auto loadgrid = [](fs::path pin,double chunk_size){
        std::ifstream in(pin);

        std::stringstream buffer;
        buffer << in.rdbuf();
        in.close();

        auto g = Grid(chunk_size);

        return g.from_xyz(buffer.str(),1/sgs::ANGSTROM);
    };

    fmt::print("generated grid, printing to file {} \n",firstout.string());
    printgrid(g,firstout);

    //main cyc

    int step = 0;

    int maxstep = 10000;

    int printstep = 200;

    auto statef = outfolder;
    statef/=("states");

    fs::create_directories(statef);


    std::unordered_map<std::shared_ptr<Atom>,kmk_data> kmk();

    auto calc_kmk_all = [](Grid& g,std::shared_ptr<grid::react::React> r){
        double mdist = r->Distance();
        for(auto iter = g.begin();iter!=g.end();) {
            auto& [vec1,atom1] = *iter.aiter;
            auto citer = g.beginFilterDistance(mdist,vec1);
            for(;!citer.Finished();) {
                auto& [vec2,atom2] = *citer.aiter;
                double chance = r->Chance(atom1,atom2,(vec2-vec1).abs());
                if(chance!=0) {
                    kmk[atom1]={r,atom1,atom2,chance};
                }
                ++citer;
            }
            ++iter;
        }
    };


    for(auto iter = g.begin();iter!=g.end();) {
        auto patom = iter.aiter->second;
        ++iter;
    }

    for(;step<maxstep;step++) {
        if(step%printstep==0) {
            auto outfile = statef;
            outfile/=(std::to_string(step)+".xyz");
            fmt::print("step {} finished, printing grid to file {} \n",step,outfile.string());
            printgrid(g,outfile);
        }



    }



}

int main() {

    grid_like_final_ex();
    return 0;
}