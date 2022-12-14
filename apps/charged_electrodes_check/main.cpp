

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
#include <chrono>



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

auto TElectrode = std::make_shared<Type>(0, __COUNTER__, "El");
auto TElectrodeR = std::make_shared<Type>(1*sgs::ELCHARGE, __COUNTER__, "Er");
auto TElectrodeL = std::make_shared<Type>(-1*sgs::ELCHARGE, __COUNTER__, "EL");
auto Oxygen = std::make_shared<Type>(0, __COUNTER__, "O");
auto Oxygen_Intersittal = std::make_shared<Type>(-2 * sgs::ELCHARGE, __COUNTER__, "OI");
auto Hafnium = std::make_shared<Type>(0, __COUNTER__, "Hf");
auto OxygenVacancy_Neutral = std::make_shared<Type>(0, __COUNTER__, "Vo");
auto OxygenVacancy_Neutral_Intersitial = std::make_shared<Type>(0, __COUNTER__, "Vi");
auto OxygenVacancy_Charged = std::make_shared<Type>(+2 * sgs::ELCHARGE, __COUNTER__, "Vp");
auto ElectrodePositive = std::make_shared<Type>(+1 * sgs::ELCHARGE, __COUNTER__, "Ep");
auto ElectrodeNegative = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "Em");

auto R1 = std::make_shared<grid::react::Puasson>(Oxygen, OxygenVacancy_Neutral_Intersitial, OxygenVacancy_Charged, Oxygen_Intersittal, 4 * sgs::ANGSTROM, sgs::ELVOLT * 3, 1e+13);
auto R2 = std::make_shared<grid::react::Puasson>(OxygenVacancy_Neutral, Oxygen_Intersittal, Oxygen_Intersittal, OxygenVacancy_Neutral, 4 * sgs::ANGSTROM, sgs::ELVOLT * 0.8, 1e+13);
auto R3 = std::make_shared<grid::react::Puasson>(OxygenVacancy_Charged, Oxygen_Intersittal, Oxygen, OxygenVacancy_Neutral_Intersitial, 4 * sgs::ANGSTROM, sgs::ELVOLT * 0.2, 1e+13);
auto R4 = std::make_shared<grid::react::Puasson>(Oxygen_Intersittal, OxygenVacancy_Neutral, OxygenVacancy_Neutral, Oxygen_Intersittal, 4 * sgs::ANGSTROM, sgs::ELVOLT * 0.8, 1e+13);

auto E1 = std::make_shared<grid::react::Puasson>(OxygenVacancy_Charged, OxygenVacancy_Neutral, OxygenVacancy_Neutral, OxygenVacancy_Charged, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.3, 1e+13);
auto E2 = std::make_shared<grid::react::Puasson>(OxygenVacancy_Charged, TElectrodeL, OxygenVacancy_Neutral, TElectrodeL, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.2, 1e+13);
auto E3 = std::make_shared<grid::react::Puasson>(OxygenVacancy_Neutral, TElectrodeR, OxygenVacancy_Charged, TElectrodeR, 5 * sgs::ANGSTROM, sgs::ELVOLT * 0.2, 1e+13);

namespace fs = std::filesystem;

class grid_runner
{
    public:
    double U_Between_Electrodes;
    double chunk_size;
    Grid g;

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
        types.push_back(OxygenVacancy_Neutral_Intersitial);
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
        g.AddType(OxygenVacancy_Neutral_Intersitial);
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
        lHfO2.add({0.78284, 0.43966, 0.7821}, OxygenVacancy_Neutral_Intersitial);
        lHfO2.add({0.21066, 0.5392, 0.22844}, OxygenVacancy_Neutral_Intersitial);

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

        double size_electrode = sgs::ANGSTROM * 8;
        double size_electrode_z = sgs::ANGSTROM * 8;
        double electrode_begin_xy = sgs::ANGSTROM * 10;
        double electrode_begin_z = electrode_end;

        Vector ElCube_s(electrode_begin_xy, electrode_begin_xy, electrode_begin_z);
        Vector ElCube_e(electrode_begin_xy + size_electrode, electrode_begin_xy + size_electrode, electrode_begin_z + size_electrode_z);

        g.ClearParallelep(ElCube_s, ElCube_e);

        int cnt = g.count();

        g.AddLattice(ElCube_s, ElCube_e, lElectrode);

        struct_end = electrode2_end;
    };

    void init_reacts_E()
    {
        E1->Name("E1");
        E2->Name("E2");
        E3->Name("E3");

        reacts.push_back(E1);
        reacts.push_back(E2);
        reacts.push_back(E3);
    }

    void init_reacts_R()
    {
        R1->Name("R1");
        R2->Name("R2");
        R3->Name("R3");
        R4->Name("R4");

        reacts.push_back(R1);
        reacts.push_back(R2);
        reacts.push_back(R3);
        reacts.push_back(R4);
    }

    void init_reacts()
    {

       init_reacts_E();
       init_reacts_R();

        


        
    };

    void init_folders(string out, string periodic_out)
    {
        outfolder = fs::absolute(out);
        fs::create_directories(outfolder);
        statef = fs::absolute(periodic_out);
        fs::create_directories(statef);
        {
            auto outfile = outfolder;
            outfile /= fmt::format("current.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
        {
            auto outfile = outfolder;
            outfile /= fmt::format("counts.txt", step);
            if(fs::exists(outfile)) fs::remove(outfile);
        }
    }

    void printgrid(fs::path pout)
    {
        std::ofstream out(pout);

        out << g.to_xyz(1 / sgs::ANGSTROM);

        out.close();
    };

    void pringgrid()
    {
        auto outfile = statef;
        outfile /= (std::to_string(step) + ".xyz");
        printgrid(outfile);
    }

    // void loadgrid(fs::path pin, double chunk_size)
    // {
    //     std::ifstream in(pin);
    //
    //     std::stringstream buffer;
    //     buffer << in.rdbuf();
    //     in.close();
    //
    //     auto g = Grid(chunk_size);
    //
    //     return g.from_xyz(buffer.str(), 1 / sgs::ANGSTROM);
    // };
    //

    void printvoltage(fs::path pout)
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
        for (auto iter = g.begin(); !(iter == g.end());)
        {
            auto &[vec, a] = *iter.aiter;
            auto &[x, y, z] = vec;
            out << fmt::format("{} {} {} {}\n", x, y, z, a->Q());
            ++iter;
        }
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
        out << step << " ";
        double sum = 0;
        for (auto& a:reacts) {
            auto ai = elsum.find(a->Name());
            if(ai!=elsum.end()) {
                double dq = ai->second/dt;
                sum+=dq;
                out << dq << " ";
            }
            else
            {
                out<<"0 ";
            }
        }
        out << sum<<std::endl;
        elsum.clear();
        dt = 0;
        out.close();
    }

    void printrcnt(fs::path pout)
    {
        std::ofstream out(pout, std::ios_base::app);
        out << step << " ";
        for (auto &[name, cnt] : react_cnt)
        {
            out << name << ":" << cnt << " ";
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
    grid_runner(double _chunk_size, double uelectrodes) : chunk_size(_chunk_size), g(_chunk_size), U_Between_Electrodes(uelectrodes), types(), reacts(), struct_end(0), step(0), maxstep(500), printstep(100), recalc_step(100), statef(""), outfolder(""), kmk(), recieved_reaction(), kmk_sum(0), Zero_field(0, 0, 0), Cond_field(uelectrodes, 0, 0), EWALD(g),react_cnt(0),elsum(0),el_begin(0),dt(0){};

    //???????????? ???????? ???????????????????????????? ?? ???????????????????????????? ??????????????????, ?????????? ?? ???????? ???????? ?????????? ?????????????????? ??????????????.

    void change_electrodes(double lel_end)
    {
        double lu = 0;
        double ru = U_Between_Electrodes;

        for (auto iter = g.begin(); !iter.Finished();)
        {
            auto &[vec, atom] = *iter.aiter;
            if (atom->Material() == TElectrode)
            {
                //atom->fixU(false);
                if (vec.z < lel_end)
                {
                    atom->Material(TElectrodeL);
                    //atom->U(lu);
                }
                else
                {

                    atom->Material(TElectrodeR);
                    //atom->U(ru);
                }
                //atom->fixU(true);
            }
            ++iter;
        }
    };

    //???????????????????????? ?????? ?????????????????????? ?????? ?????????????????????????? ??????????-?????????? ?? ?????????? ??????????????
    void calc_kmk_all(std::shared_ptr<grid::react::React> r)
    {
        int counts = 0;
        double mdist = r->Distance();
        int cnt = 0;
        for (auto iter = g.begin(); !iter.Finished();)
        {
            auto &vec1 = iter.aiter->first;
            auto atom1 = iter.aiter->second;
            auto citer = g.beginFilterDistance(mdist, vec1);

            for (; !citer.Finished();)
            {

                auto &vec2 = citer.aiter->first;
                auto atom2 = citer.aiter->second;
                double chance = r->Chance(atom1, atom2, (vec2 - vec1).abs());
                if (chance > std::abs(kmk_sum*1e-16))
                {
                    counts++;
                    auto hiter = kmk.find(atom1);
                    kmk_sum += chance;
                    if (!(kmk_sum == kmk_sum))
                    {
                        int *a = 0;
                        assert(*a == 0);
                    }
                    kmk.insert({atom1, kmk_data{r, atom1, atom2, iter.aiter->first, citer.aiter->first, chance}});
                    recieved_reaction.insert({atom2, atom1});
                }
                //fmt::print("{:e} {:e} {}\n",vec1,vec2,chance);
                ++citer;
            }
            ++iter;
        }
        fmt::print("finded reactions {} :{}\n", r->Name(), counts);
    };

    //?????????????????????????? ???????? ???? ?????????? ????????????,
    //?????? ?????? ?????????????????????????? ???????? ?????????? ???????????? ?????????????? ?????????????? ??????????

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

    //?????????????????? ???????????? ?????? ?????????????????????? ?????????? ?? ??.??
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
        //?????????????????? ????????????????????
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
                    Cond_field.Apply(pos, atom);
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
                        Cond_field.Apply(pos, atom);
                        atom->U(atom->U() + EWALD.calc_a(pos));
                    }
                }
                if (delta2 != 0)
                {

                    EWALD.add_chunk_q(p2, delta2);
                    for (auto [pos, atom] : *chunk1)
                    {
                        Zero_field.Apply(pos, atom);
                        Cond_field.Apply(pos, atom);
                        atom->U(atom->U() + EWALD.calc_a(pos));
                    }
                }
            }
        }
        //????, ?????? ?????????????????????? ???????????? ???????????????????????????? ???????????? ???????????? ???? ???????? ?? ???????????? ??????????????????????.
        if (swap)
        {
            double t1 = a1->T();
            a1->T(a2->T());
            a2->T(t1);
        }

        //???????????? ?????????? ?????????????????? ?????????? ?????????? ?????? ???????????????????? ??????????????, ???? ?????? ?????????????? ???????????? ?????? ??????????????

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
        EWALD = field::ewald_hack(g);

        change_electrodes(struct_end / 2);


        Zero_field.Apply(g);

        EWALD.calc_all();

        printvoltage("voltage.txt");

};

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
    double Chunk_size = settings.GetReal("calculation","chunk_size",1);
    class grid_runner run_this_thing_please(Chunk_size, U_between_electrodes);

    run_this_thing_please.init_types();
    run_this_thing_please.init_reacts();
    string outfile = settings.Get("folders","output","./results");
    string outfile_periodic = settings.Get("folders","periodic_output","./results/periodic");
    run_this_thing_please.init_folders(outfile,outfile_periodic);

    run_this_thing_please.init_structure();

    run_this_thing_please.maxstep = settings.GetInteger("calculation","maxstep",10000);
    run_this_thing_please.recalc_step = settings.GetInteger("calculation","recalc_step",100);
    run_this_thing_please.printstep = settings.GetInteger("calculation","printstep",100);

    std::unordered_map<std::shared_ptr<Type>,int> cts;
    for(auto a = run_this_thing_please.g.begin();!a.Finished();) {
        auto atom = a.aiter->second;

        auto materail = atom->Material();
        auto pp = cts.find(materail);
        if(pp!=cts.end()) {
            pp->second++;
        }
        else
        {
            cts[materail] = 1;
        }
        ++a;
    }
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