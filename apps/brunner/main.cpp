

#include <iostream>

#include <sgs.hpp>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <functional>
#include <memory>
#include <filesystem>
#include <list>
#include <unordered_map>
#include <cstdlib>
#include <random>
#include <limits>

#ifndef M_PI
#define M_PI 3.1415
#endif

#include <fmt/format.h>

#include <geometry/vector.hpp>
#include <geometry/geometry.hpp>
#include <grid/cubicchunk.hpp>
#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/lattice.hpp>
#include <grid/react/react.hpp>
#include <grid/react/standart.hpp>
#include <grid/grid.hpp>
#include <sgs.hpp>
#include <math/modulo.hpp>
#include <field/ewald_hack.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>
#include <chrono>

#include <assertions.h>


#include "IniReader.hpp"

using geometry::Geometry;
using geometry::Vector;
using grid::Grid;
using grid::Lattice;
using grid::atom::Atom;
using grid::atom::Type;
using std::string;

struct kmk_data
{
    std::shared_ptr<grid::react::React> react;
    std::shared_ptr<grid::atom::Atom> f;
    std::shared_ptr<grid::atom::Atom> s;
    geometry::Vector fp;
    geometry::Vector sp;
    double chance;
};

auto TElectrode = std::make_shared<Type>(0, __COUNTER__, "El",0);
auto TElectrodeR = std::make_shared<Type>(0, __COUNTER__, "Elr",0);
auto TElectrodeL = std::make_shared<Type>(0, __COUNTER__, "Ell",0);
auto Oxygen = std::make_shared<Type>(0, __COUNTER__, "O",6.4*powf(sgs::ANGSTROM,3));
auto Oxygen_Intersittal = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "OI",6.4*powf(sgs::ANGSTROM,3));
auto Hafnium = std::make_shared<Type>(0, __COUNTER__, "Hf",21.88*powf(sgs::ANGSTROM,3));
auto OxygenVacancy_Neutral = std::make_shared<Type>(0, __COUNTER__, "Vo",6.4*powf(sgs::ANGSTROM,3));
auto IntersitialPosition = std::make_shared<Type>(0, __COUNTER__, "Ip",6.4*powf(sgs::ANGSTROM,3));
auto OxygenVacancy_Charged = std::make_shared<Type>(+1 * sgs::ELCHARGE, __COUNTER__, "Vp",6.4*powf(sgs::ANGSTROM,3));
auto ElectrodePositive = std::make_shared<Type>(+1 * sgs::ELCHARGE, __COUNTER__, "Elp",0);
auto ElectrodeNegative = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "Eln",0);

auto R1 = std::make_shared<grid::react::Standart>(Oxygen, IntersitialPosition, OxygenVacancy_Charged, Oxygen_Intersittal, 5 * sgs::ANGSTROM, sgs::ELVOLT * 7, 1e+13);
auto R2 = std::make_shared<grid::react::Standart>(OxygenVacancy_Neutral, Oxygen, Oxygen, OxygenVacancy_Neutral, 5 * sgs::ANGSTROM, sgs::ELVOLT * 4, 1e+13);
auto R22 = std::make_shared<grid::react::Standart>(OxygenVacancy_Charged, Oxygen, Oxygen, OxygenVacancy_Charged, 5 * sgs::ANGSTROM, sgs::ELVOLT * 4, 1e+13);
auto R3 = std::make_shared<grid::react::Standart>(OxygenVacancy_Charged, Oxygen_Intersittal, Oxygen, IntersitialPosition, 5 * sgs::ANGSTROM, sgs::ELVOLT * 1.13, 1e+13);
auto R4 = std::make_shared<grid::react::Standart>(Oxygen_Intersittal, IntersitialPosition, IntersitialPosition, Oxygen_Intersittal, 5 * sgs::ANGSTROM, sgs::ELVOLT * 1, 1e+13);

auto E1 = std::make_shared<grid::react::Standart>(OxygenVacancy_Charged, OxygenVacancy_Neutral, OxygenVacancy_Neutral, OxygenVacancy_Charged, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E2 = std::make_shared<grid::react::Standart>(OxygenVacancy_Charged, TElectrodeL, OxygenVacancy_Neutral, TElectrodeL, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E3 = std::make_shared<grid::react::Standart>(OxygenVacancy_Neutral, TElectrodeR, OxygenVacancy_Charged, TElectrodeR, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E4 = std::make_shared<grid::react::Standart>(OxygenVacancy_Charged, TElectrodeR, OxygenVacancy_Neutral, TElectrodeR, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E5 = std::make_shared<grid::react::Standart>(OxygenVacancy_Neutral, TElectrodeL, OxygenVacancy_Charged, TElectrodeL, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);

double ATOM_E = 0;
namespace fs = std::filesystem;

class grid_runner
{
    public:
    double U_Between_Electrodes;
    double chunk_size;
    Grid g;
    INIReader settings;

    std::list<std::shared_ptr<Type>> types;
    std::list<std::shared_ptr<grid::react::React>> reacts;

    // TODO move this into another class
    std::unordered_multimap<std::shared_ptr<Atom>, kmk_data> kmk;
    std::unordered_multimap<std::shared_ptr<Atom>, std::shared_ptr<Atom>> recieved_reaction;
    std::unordered_map<std::string, double> elsum;
    std::unordered_map<std::string, int> react_cnt;

    double dt;

    double struct_end;

    fs::path statef;
    fs::path outfolder;

    // TOOLS

    field::Equal Zero_field;
    field::ZCondenser Cond_field;
    field::ewald_hack EWALD;

    double kmk_sum;
    double el_begin;

    size_t step;

public:
    size_t maxstep;
    size_t printstep;
    size_t recalc_step;

    void init_types()
    {

        types.push_back(TElectrode);
        types.push_back(TElectrodeL);
        types.push_back(TElectrodeR);
        types.push_back(Oxygen);
        types.push_back(Oxygen_Intersittal);
        types.push_back(Hafnium);
        types.push_back(OxygenVacancy_Neutral);
        types.push_back(IntersitialPosition);
        types.push_back(OxygenVacancy_Charged);
        types.push_back(ElectrodePositive);
        types.push_back(ElectrodeNegative);

        g.AddType(TElectrode);
        g.AddType(TElectrodeR);
        g.AddType(TElectrodeL);
        g.AddType(Oxygen);
        g.AddType(Oxygen_Intersittal);
        g.AddType(Hafnium);
        g.AddType(OxygenVacancy_Neutral);
        g.AddType(IntersitialPosition);
        g.AddType(OxygenVacancy_Charged);
        g.AddType(ElectrodePositive);
        g.AddType(ElectrodeNegative);
    }

    void init_structure()
    {
        double ela = 2.866 * sgs::ANGSTROM;

        Geometry geoElectrode(Vector(ela, 0, 0), Vector(0, ela, 0), Vector(0, 0, ela));

        Lattice lElectrode(geoElectrode);
        lElectrode.add(Vector(0, 0, 0), TElectrode);


        Lattice lElectrodeR(geoElectrode);
        lElectrodeR.add(Vector(0, 0, 0), TElectrodeR);

        Lattice lElectrodeL(geoElectrode);
        lElectrodeL.add(Vector(0, 0, 0), TElectrodeL);

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
        lHfO2.add({0.78284, 0.43966, 0.7821}, IntersitialPosition);
        lHfO2.add({0.21066, 0.5392, 0.22844}, IntersitialPosition);

        double size_x = 30 * sgs::ANGSTROM;
        double size_y = 30 * sgs::ANGSTROM;

        double dist_electrode = 0.3 * sgs::ANGSTROM;

        double electrode_end = 1 * sgs::ANGSTROM;

        double oxyde_begin = electrode_end + dist_electrode;

        double oxyde_end = 14 * sgs::ANGSTROM + oxyde_begin;


        el_begin = oxyde_end+dist_electrode;


        double electrode2_end = 2 * sgs::ANGSTROM + el_begin;

        g.AddLattice({0., 0, 0}, Vector(size_x, size_y, electrode_end), lElectrode);
        g.AddLattice({0, 0, oxyde_begin}, Vector(size_x, size_y, oxyde_end), lHfO2);
        g.AddLattice({0, 0, el_begin}, Vector(size_x, size_y, electrode2_end), lElectrode);

        double size_electrode = sgs::ANGSTROM * 5;
        double size_electrode_z = sgs::ANGSTROM * 5;
        double electrode_begin_xy = sgs::ANGSTROM * 10;
        double electrode_begin_z = electrode_end;

        Vector ElCube_s(electrode_begin_xy, electrode_begin_xy, electrode_begin_z);
        Vector ElCube_e(electrode_begin_xy + size_electrode, electrode_begin_xy + size_electrode, electrode_begin_z + size_electrode_z);

        g.ClearParallelep(ElCube_s, ElCube_e);

        int cnt = g.count();

        g.AddLattice(ElCube_s, ElCube_e, lElectrode);

        struct_end = electrode2_end;
    };

    void loadstructure(std::string file) {
        loadgrid(file,chunk_size);
    }

    void init_reacts_E()
    {
        E1->Name("E1");
        E2->Name("E2");
        E3->Name("E3");
        E4->Name("E4");
        E5->Name("E5");

        reacts.push_back(E1);
        reacts.push_back(E2);
        reacts.push_back(E3);
        reacts.push_back(E4);
        reacts.push_back(E5);
    }

    void init_reacts_R()
    {
        R1->Name("R1");
        R2->Name("R2");
        R22->Name("R22");
        R3->Name("R3");
        R4->Name("R4");

        reacts.push_back(R1);
        reacts.push_back(R2);
        reacts.push_back(R22);
        reacts.push_back(R3);
        reacts.push_back(R4);
    }

    void init_reacts()
    {

       init_reacts_R();
       init_reacts_E();

        


        
    };

    void init_folders(string out, string periodic_out)
    {
        outfolder = fs::absolute(out);
        statef = fs::absolute(periodic_out);
        {
            auto outfile = outfolder;
            outfile /= fmt::format("current.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
        {
            auto outfile = outfolder;
            outfile /= fmt::format("current_reacts.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
        {
            auto outfile = outfolder;
            outfile /= fmt::format("counts.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
        if(settings.GetBoolean("folder","removeout",false)) {
            if(fs::exists(outfolder)) fs::remove_all(statef);
        }
        fs::create_directories(outfolder);
        fs::create_directories(statef);
    }

    void printgrid(fs::path pout)
    {
        std::ofstream out(pout);

        out << g.to_xyz(1 / sgs::ANGSTROM);

        out.close();
    };

    auto printgrid_simple(fs::path prefix="") 
    {
        auto outfile = statef;
        outfile /= (prefix.string()+(std::to_string(step) + ".xyz"));
        fmt::print("step {} finished, printing grid to file {} \n", step, outfile.string());
        printgrid(outfile);
        return outfile;
    }

    void pringgrid()
    {
        auto outfile = statef;
        outfile /= (std::to_string(step) + ".xyz");
        printgrid(outfile);
    }

    void loadgrid(std::string& path, double chunk_size)
    {
        std::ifstream in(path);
    
        return g.from_xyz(in, 1 / sgs::ANGSTROM);
    };
    

    void printvoltage(fs::path pout)
    {
        std::ofstream out(pout);
        out << "x y z c\n";
        g.for_each([&](auto& pos,auto atom) mutable
        {
            auto &[x, y, z] = pos;
            out << fmt::format("{} {} {} {}\n", x, y, z, atom->U());
        });
        out.close();
    };
    void printvoltage()
    {

        auto outfile = statef;
        outfile /= fmt::format("voltage_{}.txt", step);
        printvoltage(outfile);
    }

    void printcharges(fs::path pout)
    {
        std::ofstream out(pout);
        out << "x y z c\n";
        g.for_each([&](auto& pos,auto atom) mutable
        {
            auto &[x, y, z] = pos;
            out << fmt::format("{} {} {} {}\n", x, y, z, atom->Q());
        });
        out.close();
    };

    void printcharges()
    {
        auto outfile = statef;
        outfile /= fmt::format("charge_{}.txt", step);
        printcharges(outfile);
    }

    void printcurrent(fs::path pout)
    {
        std::ofstream out(pout, std::ios_base::app);
        out << step;
        out << std::setprecision(6);
        double sum = 0;
        
        for (auto& a:reacts) {
            auto ai = elsum.find(a->Name());
            double dq = 0;
            if(ai!=elsum.end()) {
                dq = ai->second/dt;
                sum+=dq;
            }
            out<<"\t"<<std::setw(11) << dq;
        }
        out << "\t"<<std::setw(11)<<sum<<std::endl;
        out.close();
    }

    void printcurrent_reacts(fs::path pout)
    {
        std::ofstream out(pout, std::ios_base::app);
        double factorc=1e-20;
        out << step;
        double sum = 0;
        for (auto& a:reacts) {
            auto ai = elsum.find(a->Name());
            if(ai!=elsum.end()) {
                double dq = ai->second/dt;
                sum+=dq;
            }
        }
        double area = g.Sizes().x*g.Sizes().y;
        sum*=(factorc/area);
        double pf = calc_PF();
        double shottky = calc_Shottky();
        double direct = calc_Direct();
        double fn = calc_FN();
        double sum_all = sum+pf+shottky+direct+fn;
        out<<std::setprecision(6);
        for (auto& current:{pf,shottky,fn,direct,sum,sum_all}) {
            out<<"\t"<<std::setw(10)<<current;
        }
        out<<std::endl;
        out.close();
    }

    void printcurrent_reacts() {
        auto outfile = outfolder;
        outfile /= fmt::format("current_reacts.txt", step);
        printcurrent_reacts(outfile);
    }

    void printrcnt(fs::path pout)
    {
        std::ofstream out(pout, std::ios_base::app);
        out << fmt::format("{:6d}",step);
        for (const auto &react : reacts)
        {
            int cnt = 0;
            auto& name = react->Name();

            const auto& iterator = react_cnt.find(name);
            if(iterator!=react_cnt.end()) {
                cnt = iterator->second;
            }
            out << " " << fmt::format("{:6d}",cnt);
        }
        out << std::endl;
        react_cnt.clear();
        out.close();
    }

    void printcurrent()
    {
        auto outfile = outfolder;
        outfile /= fmt::format("current.txt", step);
        printcurrent(outfile);
    }

    void printrcnt()
    {
        auto outfile = outfolder;
        outfile /= fmt::format("counts.txt", step);
        printrcnt(outfile);
    }
    grid_runner(double _chunk_size, double uelectrodes,INIReader _settings) : chunk_size(_chunk_size), g(_chunk_size), U_Between_Electrodes(uelectrodes), types(), reacts(), struct_end(0), step(0), maxstep(500), printstep(100), recalc_step(100), statef(""), outfolder(""), kmk(), recieved_reaction(), kmk_sum(0), Zero_field(0, 0, 0), Cond_field(uelectrodes, 0, 0), EWALD(g),react_cnt(0),elsum(0),el_begin(0),dt(0),settings(_settings){};

    //меняет типы положительному и отрицательному электроду, чтобы с ними было проще указывать реакции.

    void change_electrodes(double lel_end)
    {
        size_t cnt_left = 0;
        size_t cnt_right = 0;
        g.for_each([&](auto& pos,auto atom) mutable
        {
            if (atom->Material() == TElectrode)
            {
                if (pos.z < lel_end)
                {
                    cnt_left++;
                }
                else
                {
                    cnt_right++;
                }
            }
        });
        auto size_vector = g.Rlim()-g.Llim();
        auto area = size_vector.x*size_vector.y;
        double dist_between_electrodes = size_vector.z;
        TElectrodeL->Q(area*U_Between_Electrodes/(cnt_left*dist_between_electrodes));
        TElectrodeR->Q(area*U_Between_Electrodes/(cnt_right*dist_between_electrodes));
        g.for_each([&](auto& pos,auto atom) mutable
        {
            if (atom->Material() == TElectrode)
            {
                if (pos.z < lel_end)
                {
                    atom->Material(TElectrodeL);
                }
                else
                {
                    atom->Material(TElectrodeR);
                }
            }

        });
    };

    //просчитывает все вероятности для кинетического монте-карло с одной реакции
    void calc_kmk_all(std::shared_ptr<grid::react::React> r)
    {
        int counts = 0;
        double mdist = r->Distance();
        int cnt = 0;
        g.for_each([&](const auto& pos1,auto atom1) mutable
        {
            g.for_each(pos1,mdist,[&](const auto& pos2,auto atom2) mutable
            {
                double chance = r->Chance(atom1, atom2, (pos2 - pos1).abs());
                if (chance > std::abs(kmk_sum*1e-16))
                {
                    counts++;
                    auto hiter = kmk.find(atom1);
                    kmk_sum += chance;
                    assert_simple(kmk_sum);
                    kmk.insert({atom1, kmk_data{r, atom1, atom2, pos1, pos2, chance}});
                    recieved_reaction.insert({atom2, atom1});
                }
            });
        });

        fmt::print("reactions found {} :{}\n", r->Name(), counts);
    };

    Vector getF(const Vector& pos){
        Vector f{0,0,0};
        g.for_each(pos,2*sgs::ANGSTROM,[&](const auto& pA, const auto& atom) mutable{
            Vector delta = pA-pos;

            const double cutoff = 1e-2*sgs::ANGSTROM;

            Vector deltam = {
                delta.x>cutoff?1/delta.x:0,
                delta.y>cutoff?1/delta.y:0,
                delta.z>cutoff?1/delta.z:0,
            };

            f+=deltam*atom->U();

        });

        return f;
    }

    double calc_Shottky() {
        double j = 0;
        constexpr double E_bandgap = sgs::ELVOLT*3.5;
        double E_r = 0.35*sgs::ELVOLT;
        double e_o = 0;
        double e_opt = 4;
        auto sizes = g.Sizes();
        double volume = sizes.x*sizes.y*sizes.z;
        g.for_each([&](const auto& pos,const auto& atom) mutable{
            double prefactor = sgs::ELCHARGE*sgs::ELMASS*powf(sgs::BOLZMAN,2)*atom->T()/(sgs::PLANCK/(4*M_PI));

            double field_prot = getF(pos).z;
            
            const double cutoff = 1e-5*sgs::VOLT/sgs::METER;

            if(field_prot>cutoff) {
                double sign = fabs(field_prot)/field_prot;
                j+=sign*prefactor*exp(-1/(sgs::BOLZMAN*atom->T())*(E_r-sqrt(powf(sgs::ELCHARGE,3)*fabs(field_prot))/(4*M_PI*e_opt)))*atom->V()/volume;
            }

        });

        
        return j;
    };

    double calc_PF() {
        double j = 0;
        constexpr double E_bandgap = sgs::ELVOLT*3.5;
        double E_f = 3.5*sgs::ELVOLT;
        double e_opt = 4;
        double mu = 0.1*sgs::METER/sgs::VOLT;
        double cf = 3/(Hafnium->V()+2*Oxygen->V());
        auto sizes = g.Sizes();
        double volume = sizes.x*sizes.y*sizes.z;
        g.for_each([&](const auto& pos,const auto& atom){
            double prefactor = sgs::ELCHARGE*mu*cf;

            double field_prot = getF(pos).z;

            const double cutoff = 1e-5*sgs::VOLT/sgs::METER;

            if(field_prot>cutoff) {
                double sign = fabs(field_prot)/field_prot;

                double usqrt = powf(sgs::ELCHARGE,3)*fabs(field_prot)/(4*M_PI*e_opt);

                double a = sgs::BOLZMAN*atom->T();
                double b = sqrt(usqrt);
                double c = E_f-b;
                double d = c/a;

                double uexp = -d;

                j+=sign*prefactor*exp(uexp)*atom->V()/volume;
            


            }

        });
        return j;
    };

    double calc_FN() {
        double j = 0;
        constexpr double E_bandgap = sgs::ELVOLT*3.5;
        double E_f = 3.5*sgs::ELVOLT;
        double e_opt = 4;
        double mu = 0.1*sgs::METER/sgs::VOLT;
        double sub_mass = sgs::ELMASS;
        double E_b = 1.2*sgs::ELVOLT;
        auto sizes = g.Sizes();
        double volume = sizes.x*sizes.y*sizes.z;
        g.for_each([&](const auto& pos,const auto& atom){

            double field_prot = getF(pos).z;
            double sign = fabs(field_prot)/field_prot;

            double prefactor = powf(sgs::ELCHARGE,3)*powf(fabs(field_prot),2)/(8*M_PI*sgs::PLANCK*E_b);
            const double cutoff = 1e-5*sgs::VOLT/sgs::METER;

            if(field_prot>cutoff) {
                j+=sign*prefactor*exp(-4*sqrt(sub_mass*2*powf(E_b,3)*2*M_PI/(3*sgs::PLANCK*sgs::ELCHARGE*fabs(field_prot))))*atom->V()/volume;
            }
        });
        return j;
    };

    double calc_Direct() {
        return 0; //TODO
    };

    //пересчитывает поле по всему объему,
    //так как пересчитывать поле после каждой реакции слишком долго

    void recalc_all_reactions()
    {
        kmk.clear();
        recieved_reaction.clear();

        kmk.reserve(60000);
        recieved_reaction.reserve(60000);
        kmk.max_load_factor(0.25);
        recieved_reaction.max_load_factor(0.25);

        kmk_sum = 0;

        for (auto &r : reacts)
        {
            calc_kmk_all(r);
        }
    };

    //описывает логику для перемещения атома и т.д
    double ChangeAtoms(const Vector &p1, const Vector &p2, std::shared_ptr<Type> t1, std::shared_ptr<Type> t2, bool swap)
    {
        double dq;
        auto a1 = g.get(p1);
        auto a2 = g.get(p2);
        if (a1 == nullptr || a2 == nullptr)
            return 0;
        double q1 = a1->Q();
        double q2 = a2->Q();
        a1->Material(t1);
        a2->Material(t2);

        double delta1 = a1->Q() - q1;
        double delta2 = a2->Q() - q2;
        dq = (delta1*(p2 - p1).z + delta2*(p1 - p2).z) / el_begin;
        //обновляем напряжения
        if (delta1 == 0 && delta2 == 0)
        {
        }
        else
        {

            auto chunk1 = g.getChunk(p1);
            auto chunk2 = g.getChunk(p2);

            if (chunk1 == chunk2)
            {

                double delta = delta1 + delta2;

                EWALD.add_chunk_q(p1, delta);
                for (auto [pos, atom] : *chunk1)
                {
                    Zero_field.Apply(pos, atom);
                    //Cond_field.Apply(pos, atom);
                    atom->U(atom->U() + EWALD.calc_a(pos));
                }
            }
            else
            {
                if (delta1 != 0)
                {

                    EWALD.add_chunk_q(p1, delta1);
                    for (auto [pos, atom] : *chunk1)
                    {
                        Zero_field.Apply(pos, atom);
                        //Cond_field.Apply(pos, atom);
                        atom->U(atom->U() + EWALD.calc_a(pos));
                    }
                }
                if (delta2 != 0)
                {

                    EWALD.add_chunk_q(p2, delta2);
                    for (auto [pos, atom] : *chunk1)
                    {
                        Zero_field.Apply(pos, atom);
                        //Cond_field.Apply(pos, atom);
                        atom->U(atom->U() + EWALD.calc_a(pos));
                    }
                }
            }
        }
        //да, при перемещинии атомов предпологается просто менять им типы и менять температуры (если есть температура, которую можно менять).
        if (swap)
        {
            double t1 = a1->T();
            a1->T(a2->T());
            a2->T(t1);
        }

        //теперь нужно посчитать новые шансы для затронутых реакций, ну или удалить лишние как минимум

        auto rng = kmk.equal_range(a1);
        for (auto iter = rng.first; iter != rng.second;)
        {
            auto next = iter;
            auto rngrec = recieved_reaction.equal_range(iter->second.s);
            for (auto iter2 = rngrec.first; iter2 != rngrec.second;)
            {
                auto nnext = iter2;
                nnext++;
                recieved_reaction.erase(iter2);
                iter2 = nnext;
            }

            kmk_sum -= iter->second.chance;
            next++;
            kmk.erase(iter);
            iter = next;
        }
        return dq;
    };

    void run()
    {
        Zero_field = field::Equal(0, 0, struct_end);
        Cond_field = field::ZCondenser(U_Between_Electrodes, 0, struct_end);
        EWALD = field::ewald_hack(g);

        printgrid_simple("initial_");
        change_electrodes(struct_end / 2);
        printgrid_simple("electrodes_updated_");
        static std::random_device dev;
        static std::mt19937 rng(dev());

        bool recalc = true;
        for (; step <= maxstep; step++)
        {
            

            if (step % recalc_step == 0 or recalc)
            {
                Zero_field.Apply(g);

                //Cond_field.Apply(g);

                EWALD.calc_all();

                recalc_all_reactions();

                recalc = false;

            }

            fmt::print("step: {}\n", step);
            if (step % printstep == 0)
            {
                auto outfile = printgrid_simple();
                fmt::print("step {} finished, printing grid to file {} \n", step, outfile.string());
                printvoltage();
                printcharges();
                printrcnt();
                printcurrent();
                printcurrent_reacts();
                dt = 0;
                elsum.clear();
            }

            std::uniform_real_distribution<double> dist(0, kmk_sum);
            fmt::print("step: {} sum: {}\n", step, kmk_sum);

            double rand_targ = dist(rng);
            double search_sum = 0;
            kmk_data react_info;
            size_t cccnt = 0;
            for (auto &[atom, data] : kmk)
            {
                //fmt::print("{:e} {:e}\n",data.sp,data.fp);
                cccnt++;
                search_sum += data.chance;
                if (search_sum >= rand_targ)
                {
                    react_info = data;
                    break;
                }
            }
            if(cccnt==0) {
                recalc = true;
                continue;
            }
            auto &react = react_info.react;

            assert_simple(react);

            double dq = ChangeAtoms(react_info.fp, react_info.sp, react->to1, react->to2, true);

            auto name = react->Name();

            auto iels = elsum.find(name);
            if (iels == elsum.end())
            {
                elsum[name] = dq;
            }
            else
            {
                iels->second += dq;
            }

            auto iels2 = react_cnt.find(name);
            if (iels2 == react_cnt.end())
            {
                react_cnt[name] = 1;
            }
            else
            {
                iels2->second += 1;
            }
            dt += 1/react_info.chance;
        }
    }
};

void grid_like_final_ex()
{
    INIReader settings("./settings.ini");

    if(settings.ParseError()!=0) {
        fmt::print("error when parsing settings file: settings.ini\n");
    }
    else {
        fmt::print("settings file loaded: settings.ini\n");
    }

    double U_between_electrodes = settings.GetReal("model","U_between_electrodes",0);

    double Chunk_size = settings.GetReal("calculation","chunk_size",1e-8);
    class grid_runner run_this_thing_please(Chunk_size, U_between_electrodes,settings);

    run_this_thing_please.init_types();
    run_this_thing_please.init_reacts();
    string outfile = settings.Get("folders","output","./results");
    string outfile_periodic = settings.Get("folders","periodic_output","./results/periodic");
    run_this_thing_please.init_folders(outfile,outfile_periodic);

    
    if(settings.GetBoolean("init","load",false)) {
        std::string file = settings.Get("init","loadfile","");
        assert_simple(file!="");
        run_this_thing_please.loadstructure(file);
    }
    else {
        run_this_thing_please.init_structure();
    }
    

    run_this_thing_please.maxstep = settings.GetInteger("calculation","maxstep",10000);
    run_this_thing_please.recalc_step = settings.GetInteger("calculation","recalc_step",100);
    run_this_thing_please.printstep = settings.GetInteger("calculation","printstep",100);

    std::unordered_map<std::shared_ptr<Type>,int> cts;
    run_this_thing_please.g.for_each([&](auto& pos,auto atom) mutable
    {
        auto materail = atom->Material();
        auto pp = cts.find(materail);
        if(pp!=cts.end()) {
            pp->second++;
        }
        else
        {
            cts[materail] = 1;
        }
    });
    
    for(auto& [material,cnt]:cts) {
        fmt::print("structure generated, type: {}: {}\n",material->Name(),cnt);
    }
    run_this_thing_please.run();
}

int main()
{
    grid_like_final_ex();
    return 0;
}