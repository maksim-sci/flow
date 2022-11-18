


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


void grid_addLattice() {
    Grid g(10.);

    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);

    l.add(Vector(1,1,1),std::make_shared<Type>(1,1));

    g.AddLattice(Vector(0.,0,0),Vector(30.,30.,30),l);

    fmt::print("{}: atoms inserted: {}\n",__FUNCTION__,g.count());
}

void grid_out() {
    Grid g(10.);

    auto t=std::make_shared<Type>(1,1,"O");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(t));

    auto gs = g.to_xyz();

    fmt::print("grid: \n{}\n",g.to_xyz());

    Grid g1(10.);

    g1.AddType(t);

    g1.from_xyz(gs);

    printf("good, Anakin, good!\n");

}

void grid_save_load() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(std::make_shared<Type>(1.,1.,"O")));



    auto gs = g.to_xyz();

    fmt::print("grid: \n{}\n",g.to_xyz());

    Grid g1(10.);

    g1.AddType(t);

    g1.from_xyz(gs);

    fmt::print("conts: {} {}\n",g.count(),g1.count());

    assert_eq(g.count(),g1.count());

}

void grid_save_load_big() {
    fmt::print("hello there!");
    Grid g(10.);



    auto t1 = std::make_shared<grid::atom::Type>(1.,1,"I");
    auto t2 = std::make_shared<grid::atom::Type>(1.,1,"O");
    auto a1 = std::make_shared<Atom>(t1);
    
    Vector p1 = Vector(100,100,100);
    Vector p2 = Vector(20,100,100);
    Vector p3 = Vector(20,20,20);

    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);

    l.add(Vector(1,1,1),t1);


    g.AddLattice(Vector(0.,0,0),Vector(4,4,4),l);

    assert_neq(g.count(),1);

    string string_g = g.to_xyz();

    std::ofstream out("C:/Users/sci/Desktop/modeltest/base_model.xyz");

    out<<string_g;
}

void grid_clear_parallelep_tst() {
    Grid g(3);


    Geometry geo(Vector(1,2,0)/sqrt(5),Vector(2,1.,1)/sqrt(5),Vector(1,0,2)/sqrt(5));
    Lattice l(geo);
    l.add({0,0,0},std::make_shared<Type>(0,0,"omg"));

    g.AddLattice(Vector(0.,0,0),Vector(10,10.,10),l);

    int cnt = g.count();
    std::ofstream out1("C:/Users/sci/Desktop/modeltest/base_modela.xyz");

    out1<<g.to_xyz();

    out1.close();


    std::ofstream out2("C:/Users/sci/Desktop/modeltest/base_modelb.xyz");

    out2<<g.to_xyz();

    out2.close();

    g.ClearParallelep(Vector(0.,0.,0.),Vector(7,7,7));


    std::ofstream out("C:/Users/sci/Desktop/modeltest/base_modelc.xyz");

    out<<g.to_xyz();
    //fmt::print(g.to_xyz());

    out.close();
}


void grid_like_final() {
    Grid g(sgs::ANGSTROM*10);



    auto TElectrode = std::make_shared<grid::atom::Type>(0,__COUNTER__,"El");
    auto Oxygen = std::make_shared<grid::atom::Type>(0,__COUNTER__,"O");
    auto Oxygen_Intersittal = std::make_shared<grid::atom::Type>(-2,__COUNTER__,"O");
    auto Hafnium = std::make_shared<grid::atom::Type>(0,__COUNTER__,"Hf");
    auto OxygenVacancy_Intersittal = std::make_shared<grid::atom::Type>(0,__COUNTER__,"Vo");
    auto OxygenVacancy = std::make_shared<grid::atom::Type>(+2,__COUNTER__,"Vo");

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

    lHfO2.add(Vector(0.724041, 0.457319, 0.292109), Hafnium);
    lHfO2.add(Vector(0.275959, 0.957319, 0.207891), Hafnium);
    lHfO2.add(Vector(0.275959, 0.542681, 0.707891), Hafnium);
    lHfO2.add(Vector(0.724041, 0.042681, 0.792109), Hafnium);
    lHfO2.add(Vector(0.551113, 0.742603, 0.022292), Oxygen);
    lHfO2.add(Vector(0.448887, 0.242603, 0.477708), Oxygen);
    lHfO2.add(Vector(0.448887, 0.257397, 0.977708), Oxygen);
    lHfO2.add(Vector(0.551113, 0.757397, 0.522292), Oxygen);
    lHfO2.add(Vector(0.932151, 0.330122, 0.652906), Oxygen);
    lHfO2.add(Vector(0.067849, 0.830122, 0.847094), Oxygen);
    lHfO2.add(Vector(0.067849, 0.669878, 0.347094), Oxygen);
    lHfO2.add(Vector(0.932151, 0.169878, 0.152906), Oxygen);
    lHfO2.add(Vector(0.522616, 0.110041, 0.782365), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.500000, 0.639959, 0.250000), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.833333, 0.000000, 0.500000), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.747211, 0.514587, 0.862533), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.666667, 0.000000, 0.166667), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.086123, 0.152080, 0.470800), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.283925, 0.625000, 0.048547), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.904740, 0.646880, 0.592346), OxygenVacancy_Intersittal);


    double size_x = 30*sgs::ANGSTROM;
    double size_y = 30*sgs::ANGSTROM;

    double electrode_z = 5*sgs::ANGSTROM;

    double oxyde_z = 30*sgs::ANGSTROM+electrode_z;

    double electrode2_z = 5*sgs::ANGSTROM+oxyde_z;


    g.AddLattice(Vector(0.,0,0),Vector(size_x,size_y,electrode_z),lElectrode);
    g.AddLattice(Vector(0,0,electrode_z),Vector(size_x,size_y,oxyde_z),lHfO2);
    g.AddLattice(Vector(0,0,oxyde_z),Vector(size_x,size_y,electrode2_z),lElectrode);

    double size_electrode = sgs::ANGSTROM*10;
    double electrode_begin_xy = sgs::ANGSTROM*7;
    double electrode_begin_z = electrode_z;

    Vector ElCube_s(electrode_begin_xy,electrode_begin_xy,electrode_begin_z);
    Vector ElCube_e(electrode_begin_xy+size_electrode,electrode_begin_xy+size_electrode,electrode_begin_z+size_electrode);

    g.ClearParallelep(ElCube_s,ElCube_e);

    g.AddLattice(ElCube_s,ElCube_e,lElectrode);

    fmt::print("adding electrode in zone: {:e} {:e}\n",ElCube_s,ElCube_e);

    string string_g = g.to_xyz(1/sgs::ANGSTROM);

    string f1 = "C:/Users/sci/Desktop/modeltest/model1.xyz";

    std::ofstream out(f1);

    fmt::print("sizes: {:e} {:e}\ncounts: {}\n",g.Llim(),g.Rlim(),g.count());
    out<<string_g;

    out.close();

    string f2 = "C:/Users/sci/Desktop/modeltest/model_loaded.xyz";
    std::ifstream in(f1);


    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();

    Grid g1(10*sgs::ANGSTROM);

    g1.AddType(TElectrode);
    g1.AddType(Oxygen);
    g1.AddType(Hafnium);
    g1.AddType(Oxygen_Intersittal);
    g1.AddType(OxygenVacancy_Intersittal);
    g1.AddType(OxygenVacancy);

    string s = buffer.str();

    g1.from_xyz(s,1/sgs::ANGSTROM);

    fmt::print("loaded new grid: {}",g1.count());

    out.open(f2,std::ios_base::out);
    out<<g1.to_xyz(1/sgs::ANGSTROM);
    out.close();



    

}



void grid_radius_iterator_cyclic() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    if(g.insert(Vector(1,1.,0.05),std::make_shared<Atom>(std::make_shared<Type>(1.,1.,"O")))==grid::InsertionResults::OK);

    g.setPeriod(Vector(1.,1.,1));
    g.Cyclic<'x'>(true);
    g.Cyclic<'y'>(true);

    

    Vector v(0,0,0);
    auto iter = g.beginFilterDistance(0.1,v);
    size_t cnt = 0;
    while(!iter.Finished()) {
        cnt++;
        ++iter;
    }
    assert_eq(cnt,1);


}

void grid_radius_iterator_ex() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    size_t b_cnt = 0;
    for (double d = 1; d<100; d+=1) {
        if(g.insert(Vector(d,1.,1.),std::make_shared<Atom>(std::make_shared<Type>(1.,1.,"O")))==grid::InsertionResults::OK) b_cnt++;
    }


    Vector v(0,0,0);
    auto iter = g.beginFilterDistance(150,v);
    size_t cnt = 0;
    while(!iter.Finished()) {
        cnt++;
        ++iter;
    }
    assert_eq(cnt,b_cnt);


}



void grid_radius_iter_filtered() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(std::make_shared<Type>(1.,1.,"O")));
    g.insert(Vector(1000.,1.,1.),std::make_shared<Atom>(std::make_shared<Type>(1.,1.,"O")));


    Vector v(0,0,0);
    auto iter = g.beginFilterDistance(100,v);
    size_t cnt = 0;
    while(!iter.Finished()) {
        assert_eq(iter.aiter->first,Vector(1,1,1));
        cnt++;
        ++iter;
    }
    assert_eq(cnt,1);


}

void grid_radius_iterator() {
    Grid g(10.);
    fmt::print("1\n");
    auto t = std::make_shared<Type>(1,1,"O");
    fmt::print("1\n");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(t));

    fmt::print("1\n");

    Vector v(0,0,0);
    auto iter = g.beginFilterDistance(100,v);
    size_t cnt = 0;
    fmt::print("1\n");

    while(!iter.Finished()) {
        assert_tst((iter.aiter->second->Material())==t);
        fmt::print("2\n");

        cnt++;
        ++iter;
    }
    assert_eq(cnt,1);
    fmt::print("3\n");


}

Geometry g({5.069186, 0.000000 ,-0.864173},{0.000000 ,5.195148, 0.000000},{0.000000, 0.000000, 5.326038});

#define printex(a) fmt::print("{}: {}\n",#a,a)
void pr(Vector a) {fmt::print("{:.6f} {:.6f} {:.6f} V\n",a.x,a.y,a.z);};

Vector ov(int n, double a, double b, double c) {Vector vcs[9] = {
    Vector(0,0,0),
    Vector(0.551113, 0.742603, 0.022292),
    Vector(0.448887, 0.242603, 0.477708),
    Vector(0.448887, 0.257397, 0.977708),
    Vector(0.551113, 0.757397, 0.522292),
    Vector(0.932151, 0.330122, 0.652906),
    Vector(0.067849, 0.830122, 0.847094),
    Vector(0.067849, 0.669878, 0.347094),
    Vector(0.932151, 0.169878, 0.152906)};
    return (vcs[n] + Vector(a,b,c));
    };


using std::make_shared;
void grid_ewald_hack_simple_demo() {
    
    Grid g(5.);

    auto t1 = make_shared<Type>(1.,1,"O");

    for(int x = 0; x<10;x++)
    for(int y = 0; y<10;y++)
    {
        g.insert({(double)x,(double)y,2.5},make_shared<Atom>(t1));
    }

    g.Cyclic<'x'>(true);
    g.Cyclic<'y'>(true);


    g.setPeriod({10.,10.,0});


    field::ewald_hack hh(g);



    field::ZCondenser cond(1,0,10);

    field::Equal eq(0,0,10);

    //cond.Apply(g);

    eq.Apply(g);

    hh.calc_all();

    auto iend = g.end();
    auto iter = g.begin();
   // std::ofstream out("c:/users/sci/desktop/dataomg1");
    for(;!(iter==iend);) {
        //out<<fmt::format("{} {} {}\n",iter.aiter->first.x,iter.aiter->first.y,iter.aiter->second->U());
        ++iter;
    }
    //out.close();



}

void grid_like_final_ex() {
    Grid g(sgs::ANGSTROM*10);



    auto TElectrode = std::make_shared<grid::atom::Type>(0,__COUNTER__,"El");
    auto Oxygen = std::make_shared<grid::atom::Type>(0,__COUNTER__,"O");
    auto Oxygen_Intersittal = std::make_shared<grid::atom::Type>(-2*sgs::ELCHARGE,__COUNTER__,"O");
    auto Hafnium = std::make_shared<grid::atom::Type>(0,__COUNTER__,"Hf");
    auto OxygenVacancy_Intersittal = std::make_shared<grid::atom::Type>(0,__COUNTER__,"Vo");
    auto OxygenVacancy = std::make_shared<grid::atom::Type>(+2*sgs::ELCHARGE,__COUNTER__,"Vo");

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

    lHfO2.add(Vector(0.724041, 0.457319, 0.292109), Hafnium);
    lHfO2.add(Vector(0.275959, 0.957319, 0.207891), Hafnium);
    lHfO2.add(Vector(0.275959, 0.542681, 0.707891), Hafnium);
    lHfO2.add(Vector(0.724041, 0.042681, 0.792109), Hafnium);
    lHfO2.add(Vector(0.551113, 0.742603, 0.022292), Oxygen);
    lHfO2.add(Vector(0.448887, 0.242603, 0.477708), Oxygen);
    lHfO2.add(Vector(0.448887, 0.257397, 0.977708), Oxygen);
    lHfO2.add(Vector(0.551113, 0.757397, 0.522292), Oxygen);
    lHfO2.add(Vector(0.932151, 0.330122, 0.652906), Oxygen);
    lHfO2.add(Vector(0.067849, 0.830122, 0.847094), Oxygen);
    lHfO2.add(Vector(0.067849, 0.669878, 0.347094), Oxygen);
    lHfO2.add(Vector(0.932151, 0.169878, 0.152906), Oxygen);
    lHfO2.add(Vector(0.522616, 0.110041, 0.782365), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.500000, 0.639959, 0.250000), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.833333, 0.000000, 0.500000), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.747211, 0.514587, 0.862533), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.666667, 0.000000, 0.166667), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.086123, 0.152080, 0.470800), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.283925, 0.625000, 0.048547), OxygenVacancy_Intersittal);
    lHfO2.add(Vector(0.904740, 0.646880, 0.592346), OxygenVacancy_Intersittal);


    double size_x = 30*sgs::ANGSTROM;
    double size_y = 30*sgs::ANGSTROM;

    double electrode_z = 5*sgs::ANGSTROM;

    double oxyde_z = 30*sgs::ANGSTROM+electrode_z;

    double electrode2_z = 5*sgs::ANGSTROM+oxyde_z;


    g.AddLattice(Vector(0.,0,0),Vector(size_x,size_y,electrode_z),lElectrode);
    g.AddLattice(Vector(0,0,electrode_z),Vector(size_x,size_y,oxyde_z),lHfO2);
    g.AddLattice(Vector(0,0,oxyde_z),Vector(size_x,size_y,electrode2_z),lElectrode);

    double size_electrode = sgs::ANGSTROM*10;
    double electrode_begin_xy = sgs::ANGSTROM*7;
    double electrode_begin_z = electrode_z;

    Vector ElCube_s(electrode_begin_xy,electrode_begin_xy,electrode_begin_z);
    Vector ElCube_e(electrode_begin_xy+size_electrode,electrode_begin_xy+size_electrode,electrode_begin_z+size_electrode);

    g.ClearParallelep(ElCube_s,ElCube_e);

    g.AddLattice(ElCube_s,ElCube_e,lElectrode);

    fmt::print("adding electrode in zone: {:e} {:e}\n",ElCube_s,ElCube_e);

    string string_g = g.to_xyz(1/sgs::ANGSTROM);

    string f1 = "C:/Users/sci/Desktop/modeltest/model1.xyz";

    



    

}

int main() {

    Vector v(1.,5.,6.);
    Vector v1(1e-7,5.,6.);

    // fmt::print("vector: {}\n",v);
    // fmt::print("vector: {:e}\n",v1);

    //printex(math::modulo(10.,1.1));
    //printex(math::modulo(-10.,1.1));
    //printex(math::modulo(-1.0,1.1));
    ////grid_addLattice();
    //grid_save_load();

    //grid_save_load_big();

    //grid_like_final();

    //grid_radius_iterator_ex();

    //grid_radius_iterator_cyclic();

    //grid_radius_iterator_ex();


    //grid_radius_iter_filtered();

    Vector o1(0.551113, 0.742603, 0.022292);
    Vector o2(0.448887, 0.242603, 0.477708);
    Vector o3(0.448887, 0.257397, 0.977708);
    Vector o4(0.551113, 0.757397, 0.522292);
    Vector o5(0.932151, 0.330122, 0.652906);
    Vector o6(0.067849, 0.830122, 0.847094);
    Vector o7(0.067849, 0.669878, 0.347094);
    Vector o8(0.932151, 0.169878, 0.152906);

    // Vector i1((o4+o7+o2+o1)/4);
    // Vector i2((o2+o3+o4+o5)/4);
    // Vector i3((o8+o7+o1+o5)/4);
    // Vector i4((o1+o5+o8+o3)/4);
    // Vector i5((o1+o6+o3+o4)/4);
    // Vector i6((o2+o4+o6+o7)/4);
    // Vector i7((o1+o6+o3+o7)/4);


    
    // pr(i1);
    // pr(i2);
    // pr(i3);
    // pr(i4);
    // pr(i5);
    // pr(i6);
    // pr(i7);

    Vector Hf1(0.724041, 0.457319, 0.292109);
    Vector Hf2(0.275959, 0.957319, 0.207891);
    Vector Hf3(0.275959, 0.542681, 0.707891);
    Vector Hf4(0.724041, 0.042681, 0.792109);


    //pr((ov(1,0,0,1)+ov(2,0,1,0)+ov(4,0,0,0)+ov(3,0,1,0)+ov(6,0,0,0)+ov(6,1,1,0))/6);
    //pr((ov(4,0,0,0)+ov(7,0,0,0)+ov(2,0,0,0)+ov(1,0,0,0)+ov(3,0,0,-1)+ov(8,0,1,0))/6);
    //pr((ov(5,0,0,0)+ov(6,1,-1,0)+ov(4,0,-1,0)+ov(7,1,-1,0)+ov(8,0,0,0)+ov(2,0,0,0))/6);
    //pr((ov(6,1,0,0)+ov(1,0,0,1)+ov(3,0,0,0)+ov(8,0,0,1)+ov(5,0,0,0)+ov(4,0,0,0))/6);
    //pr((ov(4,0,0,0)+ov(1,0,0,0)+ov(3,0,1,-1)+ov(8,0,1,0)+ov(2,0,1,0)+ov(6,1,0,-1))/6);
    //pr((ov(7,1,-1,1)+ov(6,1,-1,1)+ov(8,0,0,1)+ov(5,0,0,1)+ov(2,1,0,1)+ov(7,1,0,1))/6);
    //pr((ov(7,1,-1,1)+ov(6,1,-1,1)+ov(8,0,0,1)+ov(5,0,0,1)+ov(2,1,0,1)+ov(7,1,0,1))/6);
    //pr((ov(6,0,1,0)+ov(7,0,1,1)+ov(3,0,1,0)+ov(1,0,1,1))/4);
    //pr((ov(5,0,0,0)+ov(7,1,0,0)+ov(6,1,0,0)+ov(4,0,0,0))/4);

    //grid_radius_iterator();
    //grid_ewald_hack_simple_demo();

    //grid_like_final_ex();

    grid_clear_parallelep_tst();
    return 0;
}