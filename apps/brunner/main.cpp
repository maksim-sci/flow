//Оставь надежду, всяк сюда входящий!
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
#include <grid/react/barrier.hpp>
#include <grid/react/tat.hpp>
#include <grid/react/ionic.hpp>
#include <grid/react/tat.hpp>
#include <grid/grid.hpp>
#include <sgs.hpp>
#include <math/modulo.hpp>
#include <field/ewald.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>
#include <chrono>
#include <field/summation.hpp>
#include <algo/kmk.hpp>

#include <assertions.h>


#include "IniReader.hpp"

using geometry::Geometry;
using geometry::Vector;
using grid::Grid;
using grid::Lattice;
using grid::atom::Atom;
using grid::atom::Type;
using std::string;

auto TElectrode = std::make_shared<Type>(0, __COUNTER__, "El",0);
auto TElectrodeR = std::make_shared<Type>(0, __COUNTER__, "Elr",0);
auto TElectrodeL = std::make_shared<Type>(0, __COUNTER__, "Ell",0);
auto TElectrodeP = std::make_shared<Type>(0, __COUNTER__, "Ep",0);
auto TElectrodeN = std::make_shared<Type>(0, __COUNTER__, "En",0);
auto Oxygen = std::make_shared<Type>(-1* sgs::ELCHARGE, __COUNTER__, "O",6.4*powf(sgs::ANGSTROM,3));
auto Oxygen_Intersittal = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "OI",6.4*powf(sgs::ANGSTROM,3));
auto Hafnium = std::make_shared<Type>(2* sgs::ELCHARGE, __COUNTER__, "Hf",21.88*powf(sgs::ANGSTROM,3));
auto OxygenVacancy_Neutral = std::make_shared<Type>(0, __COUNTER__, "Vo",6.4*powf(sgs::ANGSTROM,3));
auto IntersitialPosition = std::make_shared<Type>(0, __COUNTER__, "Ip",6.4*powf(sgs::ANGSTROM,3));
auto OxygenVacancy_Charged = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "Vo",6.4*powf(sgs::ANGSTROM,3));


auto R1 = std::make_shared<grid::react::ionic>(Oxygen, IntersitialPosition, OxygenVacancy_Neutral, Oxygen_Intersittal, 3 * sgs::ANGSTROM, sgs::ELVOLT * 7, 1e+13);
auto R21 = std::make_shared<grid::react::ionic>(OxygenVacancy_Neutral, Oxygen, Oxygen, OxygenVacancy_Neutral, 3 * sgs::ANGSTROM, sgs::ELVOLT * 4, 1e+13);
auto R22 = std::make_shared<grid::react::ionic>(OxygenVacancy_Charged, Oxygen, Oxygen, OxygenVacancy_Charged, 3 * sgs::ANGSTROM, sgs::ELVOLT * 4, 1e+13);
auto R3 = std::make_shared<grid::react::ionic>(OxygenVacancy_Neutral, Oxygen_Intersittal, Oxygen, IntersitialPosition, 3 * sgs::ANGSTROM, sgs::ELVOLT * 1.13, 1e+13);
auto R4 = std::make_shared<grid::react::ionic>(Oxygen_Intersittal, IntersitialPosition, IntersitialPosition, Oxygen_Intersittal, 3 * sgs::ANGSTROM, sgs::ELVOLT * 1, 1e+13);

auto E1 = std::make_shared<grid::react::ionic>(OxygenVacancy_Charged, OxygenVacancy_Neutral, OxygenVacancy_Neutral, OxygenVacancy_Charged, 3 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E21 = std::make_shared<grid::react::ionic>(OxygenVacancy_Charged, TElectrodeP, OxygenVacancy_Neutral, TElectrodeP, 3 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E31 = std::make_shared<grid::react::ionic>(OxygenVacancy_Neutral, TElectrodeN, OxygenVacancy_Charged, TElectrodeN, 3 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E22 = std::make_shared<grid::react::ionic>(OxygenVacancy_Charged, TElectrodeN, OxygenVacancy_Neutral, TElectrodeN, 3 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E32 = std::make_shared<grid::react::ionic>(OxygenVacancy_Neutral, TElectrodeP, OxygenVacancy_Charged, TElectrodeP, 3 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);


namespace fs = std::filesystem;

class basic_runner
{
    public:
    double U_Between_Electrodes;
    double chunk_size{1};
    Grid g;
    INIReader settings;

    std::list<std::shared_ptr<Type>> types;

    // TODO move this into another class
    double elsum{0};

    algo::kmk kmk_ionic;
    algo::kmk kmk_electron;

    double struct_end;

    fs::path statef{""};
    fs::path outfolder{""};

    // TOOLS

    field::Equal Zero_field;
    field::ZCondenser Cond_field;
    field::ewald EWALD;

    double elcharge_factor;

    double el_begin{0};

    size_t step{0};

public:
    size_t maxstep{0};
    size_t printstep{999999};
    size_t recalc_step{1};
    size_t calc_current{5000};

    void init_types()
    {

        types.push_back(TElectrode);
        types.push_back(TElectrodeL);
        types.push_back(TElectrodeR);
        types.push_back(TElectrodeP);
        types.push_back(TElectrodeN);
        types.push_back(Oxygen);
        types.push_back(Oxygen_Intersittal);
        types.push_back(Hafnium);
        types.push_back(OxygenVacancy_Neutral);
        types.push_back(IntersitialPosition);
        types.push_back(OxygenVacancy_Charged);

        g.AddType(TElectrode);
        g.AddType(TElectrodeR);
        g.AddType(TElectrodeL);
        g.AddType(TElectrodeP);
        g.AddType(TElectrodeN);
        g.AddType(Oxygen);
        g.AddType(Oxygen_Intersittal);
        g.AddType(Hafnium);
        g.AddType(OxygenVacancy_Neutral);
        g.AddType(IntersitialPosition);
        g.AddType(OxygenVacancy_Charged);
    }

    void init_structure()
    {
        double ela = 2.866 * sgs::ANGSTROM;

        Geometry geoElectrode(Vector(ela, 0, 0), Vector(0, ela, 0), Vector(0, 0, ela));

        Lattice lElectrode(geoElectrode);
        lElectrode.add(Vector(0, 0, 0), TElectrode);

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
        lHfO2.add({0.61202, 0.9786, 0.39605}, IntersitialPosition);
        lHfO2.add({0.38798, 0.0213, 0.60394}, IntersitialPosition);

        double size_x = settings.GetReal("model","size_x",30) * sgs::ANGSTROM;
        double size_y = settings.GetReal("model","size_y",30) * sgs::ANGSTROM;
        double size_z = settings.GetReal("model","size_z",20) * sgs::ANGSTROM;

        double dist_electrode = 0.1 * sgs::ANGSTROM;

        double electrode_end = 0.3 * sgs::ANGSTROM;

        double oxyde_begin = electrode_end + dist_electrode;

        double oxyde_end = settings.GetReal("model","oxyde_height",14)* sgs::ANGSTROM + oxyde_begin;


        el_begin = oxyde_end+dist_electrode;


        double electrode2_end = 3 * sgs::ANGSTROM + el_begin;

        g.AddLattice({0., 0, 0}, Vector(size_x, size_y, electrode_end), lElectrode);
        g.AddLattice({0, 0, oxyde_begin}, Vector(size_x, size_y, oxyde_end), lHfO2);
        g.AddLattice({0, 0, el_begin}, Vector(size_x, size_y, electrode2_end), lElectrode);

        double size_electrode = sgs::ANGSTROM * settings.GetReal("model","size_electrode",5);
        double size_electrode_z = sgs::ANGSTROM * settings.GetReal("model","size_electrode_z",5);
        double electrode_begin_xy = sgs::ANGSTROM * settings.GetReal("model","electrode_begin_xy",10);
        double electrode_begin_z = electrode_end;

        Vector ElCube_s(electrode_begin_xy, electrode_begin_xy, electrode_begin_z);
        Vector ElCube_e(electrode_begin_xy + size_electrode, electrode_begin_xy + size_electrode, electrode_begin_z + size_electrode_z);

        g.ClearParallelep(ElCube_s, ElCube_e);

        int cnt = g.count();

        g.AddLattice(ElCube_s, ElCube_e, lElectrode);

        struct_end = electrode2_end;

        g.Cyclic<'x'>(true);
        g.Cyclic<'y'>(true);
    };

    void loadstructure(std::string file) {
        loadgrid(file);

        g.for_each([&](const auto& pos, const auto& atom){
            if(
                atom->Material()==TElectrodeR||
                atom->Material()==TElectrodeL||
                atom->Material()==TElectrodeP||
                atom->Material()==TElectrodeN
            ) {
                atom->Material(TElectrode);
            }
        });
    }

    void init_reacts_E()
    {
        E1->Name("E1");
        E21->Name("E21");
        E31->Name("E31");
        E22->Name("E22");
        E32->Name("E32");
        kmk_electron.add(E1);
        kmk_electron.add(E21);
        kmk_electron.add(E31);
        kmk_electron.add(E22);
        kmk_electron.add(E32);

        // reacts.push_back(E1);
        // reacts.push_back(E21);
        // reacts.push_back(E31);
        // reacts.push_back(E22);
        // reacts.push_back(E32);
    }

    void init_reacts_R()
    {
        R1->Name("R1");
        R21->Name("R21");
        R22->Name("R22");
        R3->Name("R3");
        R4->Name("R4");

        kmk_ionic.add(R1);
        kmk_ionic.add(R21);
        kmk_ionic.add(R21);
        kmk_ionic.add(R3);
        kmk_ionic.add(R4);
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
        try{
            auto outfile = outfolder;
            outfile /= fmt::format("current.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);

        } 
        catch(std::exception e) {
            fmt::print("{} :{},'{}'\n",__FILE__,__LINE__,e.what());
        }
        try {
            auto outfile = outfolder;
            outfile /= fmt::format("current_reacts.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
        catch(std::exception e) {
            fmt::print("{} :{},'{}'\n",__FILE__,__LINE__,e.what());
        }
        try {
            auto outfile = outfolder;
            outfile /= fmt::format("counts.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
        catch(std::exception e) {
            fmt::print("{} :{},'{}'\n",__FILE__,__LINE__,e.what());
        }
        try {
            if(settings.GetBoolean("folders","removeout",false)) {
                if(fs::exists(outfolder)) fs::remove_all(statef);
            }
        }
        catch(std::exception e) {
            fmt::print("{}:{},'{}'\n",__FILE__,__LINE__,e.what());
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

    void loadgrid(std::string& path)
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
            out << fmt::format("{} {} {} {}\n", x, y, z, atom->U()/sgs::VOLT);
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
            out << fmt::format("{} {} {} {} {}\n", x, y, z, atom->Q()/sgs::ELCHARGE,atom->Material()->Name());
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
        
        out << "\t"<<std::setw(11)<<elsum/kmk_electron.Time()<<std::endl;
        out.close();
    }

    void printcurrent_reacts(fs::path pout)
    {
        std::ofstream out(pout, std::ios_base::app);
        double factorc=1e-20;
        out << step;
        double sum = elsum;
        double area = g.Sizes().x*g.Sizes().y;
        sum*=(1/area);
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


    void printcurrent()
    {
        auto outfile = outfolder;
        outfile /= fmt::format("current.txt", step);
        printcurrent(outfile);
    }

    // void printrcnt(fs::path pout)
    // {
    //     std::ofstream out(pout, std::ios_base::app);
    //     out << fmt::format("{:6d}",step);
    //     for (const auto &react : reacts)
    //     {
    //         int cnt = 0;
    //         auto& name = react->Name();

    //         const auto& iterator = react_cnt.find(name);
    //         if(iterator!=react_cnt.end()) {
    //             cnt = iterator->second;
    //         }
    //         out << " " << fmt::format("{:6d}",cnt);
    //     }
    //     out << std::endl;
    //     react_cnt.clear();
    //     out.close();
    // }

    // void printrcnt()
    // {
    //     auto outfile = outfolder;
    //     outfile /= fmt::format("counts.txt", step);
    //     printrcnt(outfile);
    // }

    basic_runner(double _chunk_size, double uelectrodes,const char* _settings) : 
    chunk_size(_chunk_size), 
    g(_chunk_size), 
    U_Between_Electrodes(uelectrodes), 
    types(), 
    kmk_ionic(&g),
    kmk_electron(&g),
    Zero_field(0, 0, 0), 
    Cond_field(uelectrodes, 0, 0), 
    settings(_settings),
    EWALD(settings.GetReal("ewald","real_cutoff",sgs::ANGSTROM*5),settings.GetReal("ewald","reciprocal_cutoff",sgs::ANGSTROM*5),settings.GetInteger("ewald","calc_size",1),settings.GetReal("ewald","sigma",1)*sgs::ANGSTROM,&g),
    elcharge_factor(settings.GetReal("model","elcharge_factor",100))
    {};

    //меняет типы положительному и отрицательному электроду, чтобы с ними было проще указывать реакции.

    void change_electrodes(double lel_end,double elcharge_factor)
    {
        double cnt = 0;
        //Специальная формула для рассчета электродов. p = d/l
        //задумана для того, чтобы имитировать скапливание заряда на дефекте
        //l = |max_electrode_defect-min_electrode_defect|+formula_epsilon
        //d = |electrode_end-pos|
        //electrode_max - максимальное расстояние от начала электрода до дефекта.
        //электроды ищутся на конце и в начале всей системы.

        constexpr double electrode_cutoff = 5*sgs::ANGSTROM; //максимальное расстояние, на котором электрод считается цельным.

        struct electrode_info {
            double begin{0};
            double end{0};
            std::shared_ptr<Type> type;
        };

        double begin_electrode_L = 0;
        double begin_electrode_R = g.Sizes().z;

        double end_electrode_L = begin_electrode_L;
        double end_electrode_R = begin_electrode_R;

        bool left_positive = U_Between_Electrodes>0;

        auto type_L = left_positive?TElectrodeP:TElectrodeN;
        auto type_R = left_positive?TElectrodeN:TElectrodeP;

        electrode_info left_electrode{begin_electrode_L,end_electrode_L,type_L};
        electrode_info right_electrode{begin_electrode_R,end_electrode_R,type_R};

        
        //шизогонический рекусивный код на лямбдах для замены просто электроды на типизированные
        const auto check_electrode_with_type = [&](const Vector& pos, electrode_info& electrode,std::shared_ptr<Type> init_type) {
            std::function<void(const Vector& pos,Grid* g)> check_impl;

            check_impl = [&electrode,&check_impl,&init_type](const Vector& pos,Grid* g)mutable {
                g->for_each(pos,electrode_cutoff,[&](const auto& pos2, const auto& atom)mutable{
                    if(atom->Material()==init_type) {
                        atom->Material(electrode.type);
                        if(std::fabs(pos2.z-electrode.begin)>std::fabs(electrode.begin-electrode.end)) {
                            electrode.end = pos2.z;
                        }
                        check_impl(pos2,g);
                    }
                });
            };
            check_impl(pos,&g);
        };

        //точки для начала поиска электродов
        Vector check_electrode_L{0,0,begin_electrode_L}; 
        Vector check_electrode_R{0,0,begin_electrode_R};

        check_electrode_with_type(check_electrode_L,left_electrode,TElectrode);
        check_electrode_with_type(check_electrode_R,right_electrode,TElectrode);

        fmt::print("left electrode: ({},{})\n",left_electrode.begin,left_electrode.end);
        fmt::print("right electrode: ({},{})\n",right_electrode.begin,right_electrode.end);

        //settings for custom formula:
        //enabled = used or not
        //epsilon: p = multiplier * (|end-pos|)/(|end-begin|+eps)

        auto check_electrode_formula = [&](std::string suffix,electrode_info& electrode,double q) {
            constexpr char section[] = "electrode_charge_custom";


            std::string section_name(section);
            section_name+=suffix;

            bool enabled = settings.GetBoolean(section_name,"enabled",false);
            double sum = 0;

            double epsilon = settings.GetReal(section_name,"eps",1*sgs::ANGSTROM);
            double multiplier = settings.GetReal(section_name,"mult",1);

            g.for_each([&](const auto& pos, auto& atom) mutable{
                if(atom->Material()==electrode.type) {
                    double q=0;
                    if(enabled)
                    {
                        q = (std::fabs(electrode.end-electrode.begin))/(epsilon+std::fabs(electrode.end-pos.z));
                    }
                    else
                    {
                        q = 1;
                    }
                    sum+=q;
                    atom->Q(q);
                }
            });

            double partial = q/sum;

            g.for_each([&](const auto& pos, auto& atom) mutable {
                if(atom->Material()==electrode.type) {
                    atom->Q(atom->Q()*partial*multiplier);
                }
            });

        };


        auto size_vector = g.Sizes();
        auto area = size_vector.x*size_vector.y;
        double dist_between_electrodes = size_vector.z;
        double charge = std::fabs(U_Between_Electrodes*area/dist_between_electrodes);
        charge*=elcharge_factor;

        if(left_positive) {
            check_electrode_formula("_P",left_electrode,charge);
            check_electrode_formula("_N",right_electrode,-charge);
        }
        else {
            check_electrode_formula("_P",right_electrode,charge);
            check_electrode_formula("_N",left_electrode,-charge);
        }





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

    //описывает логику для перемещения атома и т.д
    //swap = true - значит атомы физически переместились
    double ChangeAtoms(const Vector &p1, const Vector &p2, std::shared_ptr<Type> t1, std::shared_ptr<Type> t2, bool swap)
    {
        double dq;
        auto a1 = g.get(p1);
        auto a2 = g.get(p2);
        if (a1 == nullptr || a2 == nullptr)
            return 0;

        double delta1,delta2;
        if(a1->Material()!=t1) {
            double q1 = a1->Q();
            a1->Material(t1);
            delta1 = a1->Q() - q1; 
            EWALD.add_charge(p1,delta1);

        }
        else {
            delta1 = 0;
        }

        if(a2->Material()!=t2) {
            double q2 = a2->Q();
            a2->Material(t2);
            double delta2 = a2->Q() - q2;
            EWALD.add_charge(p2,delta2);
        }
        else {
            delta2 = 0;
        }

        dq = (delta1*(p2 - p1).z + delta2*(p1 - p2).z) / el_begin;
        if (swap)
        {
            double t1 = a1->T();
            a1->T(a2->T());
            a2->T(t1);
        }

        return dq;
    };

    void run()
    {
        Zero_field = field::Equal(0, 0, struct_end);
        Cond_field = field::ZCondenser(U_Between_Electrodes, 0, struct_end);

        printgrid_simple("initial_");
        change_electrodes(struct_end / 2,elcharge_factor);
        printgrid_simple("electrodes_updated_");

        {
            fmt::print("initializing field!\n");
            Zero_field.Apply(g);

            EWALD.cache_grid_data();
            EWALD.apply();
            //print
            auto outfile = statef;
            outfile/=fmt::format("voltage_initial.txt");
            printvoltage(outfile);
            fmt::print("field initialization completed\n");
        }

        bool recalc = true;
        bool print = true;
        bool print_current = false;
        for (; step <= maxstep;)
        {
            if(step % recalc_step == 0) recalc = true;

            if(step % printstep == 0) print = true;

            if(step%calc_current==0) print_current = true;

            if(print) {
                printcharges();
                printvoltage();
            }

            if (recalc)
            {

                printf("recalculating reactions!\n");
                kmk_ionic.recalc();

                recalc = false;

            }

            fmt::print("step: {}\n", step);
            if (print)
            {
                auto outfile = printgrid_simple();
                fmt::print("step {} finished, printing grid to file {} \n", step, outfile.string());
                printvoltage();
                elsum=0;
            }

            if(print_current) {
                printcurrent();
                printcurrent_reacts();
            }

            fmt::print("step: {} sum: {}\n", step, kmk_ionic.Sum());

            double search_sum = 0;
            
            auto [flag,react] = kmk_ionic.findAndProcessReact();
            if(not flag) {
                fmt::print("error: unable to find reaction at step {}\n",step);
            }

            double dq = ChangeAtoms(react.fp, react.sp, react.r->to1, react.r->to2, true);

            elsum+=dq;

            step++;
            print = false;
        }
        printf("run ended successfully, final step: %zu\n",step);
    }
};

void grid_like_final_ex()
{
    const char* settings_path = "./settings.ini";

    INIReader settings(settings_path);

    if(settings.ParseError()!=0) {
        fmt::print("error when parsing settings file: settings.ini\n");
    }
    else {
        fmt::print("settings file loaded: settings.ini\n");
    }

    double U_between_electrodes = -settings.GetReal("model","U_between_electrodes",0)*sgs::VOLT;

    double Chunk_size = settings.GetReal("calculation","chunk_size",1e-8);
    class basic_runner runner(Chunk_size, U_between_electrodes,settings_path);

    runner.init_types();
    runner.init_reacts();

    
    string outfile = settings.Get("folders","output","./results");
    string outfile_periodic = settings.Get("folders","periodic_output","./results/periodic");
    runner.init_folders(outfile,outfile_periodic);

    
    if(settings.GetBoolean("init","load",false)) {
        std::string file = settings.Get("init","loadfile","");
        assert_simple(file!="");
        runner.loadstructure(file);
    }
    else {
        runner.init_structure();
    }

    runner.g.Cyclic<'x'>(true);
    runner.g.Cyclic<'y'>(true);

    if(settings.GetBoolean("boundaries","used",false)) {
        double size_y = settings.GetReal("boundaries","size_y",-1)*sgs::ANGSTROM;
        double size_x = settings.GetReal("boundaries","size_x",-1)*sgs::ANGSTROM;

        auto sizes = runner.g.Sizes();
        bool change_sizes = false;

        if(size_x>0) {
            change_sizes = true;
            sizes.x = size_x;
        }
        if(size_y>0) {
            change_sizes = true;
            sizes.y = size_y;
        }

        if(change_sizes) {
            runner.g.setPeriod(sizes);
        }



    }
    

    runner.maxstep = settings.GetInteger("calculation","maxstep",10000);
    runner.recalc_step = settings.GetInteger("calculation","recalc_step",100);
    runner.calc_current = settings.GetInteger("calculation","calc_current",5000);
    runner.printstep = settings.GetInteger("calculation","printstep",100);

    std::unordered_map<std::shared_ptr<Type>,int> cts;
    runner.g.for_each([&](auto& pos,auto atom) mutable
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
    runner.run();
}

int main()
{
    grid_like_final_ex();
    return 0;
}