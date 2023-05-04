


#include <iostream>
#include <string>
#include <map>
#include <fmt/format.h>
#include <algorithm>
#include <unordered_map>

#include <sgs.hpp>
#include "macro_switch.hpp"
#include <geometry/vector.hpp>
#include <geometry/geometry.hpp>
#include <grid/cubicchunk.hpp>
#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/lattice.hpp>
#include <grid/react/react.hpp>
#include <grid/react/barrier.hpp>
#include <grid/grid.hpp>
#include <math/modulo.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>
#include <field/ewald_hack.hpp>

#include <filesystem>
#include <fstream>

#include <grid/react/barrier.hpp>

#include <assertions.h>

using std::string;
using geometry::Vector;
using geometry::Geometry;
using grid::chunk::CubicChunk;



int test_basic_sgs() {
    if(!(sgs::JOULE==sgs::ERG*1e7)) {
        return -1;
    }

    assert_eq(sgs::JOULE,sgs::ERG*1e+7);

    return 0;
}

int all_work() {
    return 0;
}

void math_modulo_test_double() {
    double d = math::modulo(10.,1.1);

    assert_tst(std::abs(d-0.1)<0.00001);
}

int test_vector_creation() {
    auto v1 = Vector();
    auto v2 = Vector(1.,1.,1.);
    auto v3 = Vector(2.,3.,4.);

    return 0;
}

int vector_add() {

    auto v1 = Vector();

    auto v2 = Vector(1.,1.,1.);
    auto v3 = Vector(2.,3.,4.);
    auto s = Vector(3.,4.,5.);

    assert_eq(v1,v1);
    assert_eq(v2+v3,s);

    return 0;
}

int vector_cmp_check() {


    auto v2 = Vector(1.,1.,1.);
    auto v3 = Vector(1.01,1.,1.);
    auto v0 = Vector(0.,0.,0.);

    assert_neq(v2,v3);
    assert_eq(v2,v2);
    assert_eq(v0,v0);

    return 0;
}

int test_must_crush() {
    throw "it is error!";
    return 0;
}

int geometry_create() {
   Geometry g(Vector(1.,0.,0.),Vector(0.,1.,0.),Vector(0.,0.,1.));
   return 0;
}

int geometry_to() {
    Geometry g1(Vector(1.,0.,0.),Vector(0.,1.,0.),Vector(0.,0.,1.));
    Geometry g2(Vector(2.,0.,0.),Vector(0.,2.,0.),Vector(0.,0.,2.));
    Geometry g3(Vector(1.,0.,0.),Vector(0.,2.,0.),Vector(0.,0.,3.));

    auto v1 = Vector(1.,2.,3.);

    if(g1.to_euclidus(v1)!=v1) return -1;
    if(g2.to_euclidus(v1)!=Vector(2.,4.,6.)) return -1;
    if(g3.to_euclidus(v1)!=Vector(1.,4.,9.)) return -1;
    


    return 0;

    
}

int geometry_from() {
    Geometry g1(Vector(1.,0.,0.),Vector(0.,1.,0.),Vector(0.,0.,1.));
    Geometry g2(Vector(2.,0.,0.),Vector(0.,2.,0.),Vector(0.,0.,2.));
    Geometry g3(Vector(1.,0.,0.),Vector(0.,2.,0.),Vector(0.,0.,3.));

    auto v1 = Vector(1.,2.,3.);
    auto v2 = Vector(2.,4.,6.);
    auto v3 = Vector(1.,4.,9.);

    assert_eq(g1.from_euclidus(v1),v1);
    assert_eq(g2.from_euclidus(v2),v1);
    assert_eq(g3.from_euclidus(v3),v1);
    


    return 0;

    
}

int vector_hash() {
    auto v1 = Vector(1.,2.,3.);
    auto v2 = Vector(2.,4.,6.);

    assert_eq(std::hash<Vector>()(v1),std::hash<Vector>()(v1));
    assert_neq(std::hash<Vector>()(v1),std::hash<Vector>()(v2));
    return 0;
}

int vector_cmp() {
    auto v1 = Vector(1.,2.,3.);
    auto v2 = Vector(2.,4.,6.);

    assert_tst(v2>v1);
    assert_tst(!(v1>v1));

    return 0;
}

int cubic_chunk_creation() {
    CubicChunk a(0,0,0,1.);
    CubicChunk b;
    CubicChunk c(Vector(0.,0.,0.),128);

    return 0;
}
using grid::atom::Atom;
using grid::atom::Type;
int cubic_chunk_insert() {
    CubicChunk a(Vector(0.,0.,0.),128);

    auto ptr = std::make_shared<Atom>(std::make_shared<Type>(1,1,""));
    auto pos = Vector(0.,0.,0.);

    a.insert(pos,ptr);

    return 0;
}



int cubic_chunk_multy() {
    CubicChunk a(Vector(0.,0.,0.),128);
    auto t = std::make_shared<grid::atom::Type>(1,1,"");

    auto aa = std::make_shared<Atom>(t);


    assert_eq(a.insert({0.,0.,0.},aa),grid::InsertionResults::OK);
    assert_eq(a.insert({.5,0.,0.},aa),grid::InsertionResults::OK);
    assert_eq(a.insert({0,.5,0.},aa),grid::InsertionResults::OK);
    assert_eq(a.insert({0,.5,0.},aa),grid::InsertionResults::RepeatedInsertion);

    return 0;
}

int cubic_chunk_insert_get() {
    CubicChunk a(Vector(0.,0.,0.),128);

    auto a1 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a2 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a3 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));


    assert_eq(a.insert(Vector(0.,0.,0.),a1),grid::InsertionResults::OK);
    assert_eq(a.insert(Vector(0.5,0.,0.),a2),grid::InsertionResults::OK);
    assert_eq(a.insert(Vector(0.,0.5,0.),a3),grid::InsertionResults::OK);

    assert_tst(a.get(Vector(0.,0.,0.))==a1);
    assert_tst(a.get(Vector(0.5,0.,0.))==a2);
    assert_tst(a.get(Vector(0,0.5,0.))==a3);

    return 0;
}


int cubic_chunk_iterate() {
    CubicChunk a(Vector(0.,0.,0.),128);

    auto a1 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a2 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a3 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));


    a.insert(Vector(0.,0.,0.),a1);
    a.insert(Vector(0.5,0.,0.),a2);
    a.insert(Vector(0.,0.5,0.),a3);


    size_t cnt = 0;
    for (auto group:a) {
        cnt++;
    }

    assert_eq(cnt,3);

    return 0;
}

void cubic_chunk_count() {
    CubicChunk a(Vector(0.,0.,0.),128);

    auto a1 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a2 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a3 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));


    a.insert(Vector(0.,0.,0.),a1);
    a.insert(Vector(0.5,0.,0.),a2);
    a.insert(Vector(0.,0.5,0.),a3);

    assert_eq(a.count(),3);
}

void cubic_chunk_erase() {
    CubicChunk a(Vector(0.,0.,0.),128);

    auto a1 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a2 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));
    auto a3 = std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1,1,""));


    a.insert(Vector(0.,0.,0.),a1);
    a.insert(Vector(0.5,0.,0.),a2);
    a.insert(Vector(0.,0.5,0.),a3);

    a.erase(Vector(0.5,0.,0.));

    assert_eq(a.count(),2);
}

using grid::atom::Type;
int atom_type_creation_and_cmp() {
    auto a = std::make_shared<grid::atom::Type>(1,1,"OXYGEN");
    auto a1 = std::make_shared<grid::atom::Type>(0,1,"ELECTRODE");
    auto b = std::make_shared<grid::atom::Type>(3);

    assert_eq(a->Name(),"OXYGEN");
    assert_eq(a1->Name(),"ELECTRODE")

    assert_neq(*a,*b);
    assert_eq(*a,*a);
    assert_eq(*a,*a1);
    return 0;
}
using grid::Lattice;
void lattice_creation_check() {


    Geometry g(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1));
    Lattice l(g);

    auto t1 = std::make_shared<grid::atom::Type>(1.,4,"ELECTRODE");
    auto t2 = std::make_shared<grid::atom::Type>(1.,5,"OXYGEN");
    auto t3 = std::make_shared<grid::atom::Type>(1.,6,"UNKNOWN");


    l.add(Vector(1.,1.,1.),t1);
    l.add(Vector(0.,1.,1.),t2);
    l.add(Vector(1.,0.,1.),t3);

    size_t cnt = 0;
    for (auto& a:l) {
        cnt++;
    }

    assert_eq(cnt,3);
}

void lattice_creation_check_incorrect_add() {
    Geometry g(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1));
    Lattice l(g);

    auto t2 = std::make_shared<grid::atom::Type>(1.,7,"OXYGEN");


    Vector vec(0.,0.,0.);


    try {
        l.add(vec,t2);
        l.add(vec,t2);
    }
    catch(string& s) {
        return;
    }
    exit(-1);
}

using grid::atom::Atom;

typedef std::shared_ptr<Atom> pAtom_t;
using std::make_shared;
void check_atom_create() {
    auto t = make_shared<Type>(0.,1,"OXYGEN");

    grid::atom::Atom a(t);
}

using grid::react::react;

void check_reaction_create() {
    auto t1=make_shared<Type>(0.,1,"OXYGEN");
    auto t2=make_shared<Type>(0.,2,"ELECTRODE");

    react r(t1,t2,t1,t2,100);
}

void check_reaction_check_atoms() {
    auto t1=make_shared<Type>(0.,1,"OXYGEN");
    auto t2=make_shared<Type>(0.,2,"ELECTRODE");


    auto a1=make_shared<Atom>(t1);

    auto a2=make_shared<Atom>(t2);

    react r(t1,t2,t1,t2,100);

    assert_tst(r.AreAtomsOk(a1,a2));
    assert_tst(!r.AreAtomsOk(a1,a1));


}
using grid::Grid;
void create_grid() {
    Grid g(128.);
}


void grid_add_Reaction() {
    Grid g(128.);

    auto t1 = std::make_shared<grid::atom::Type>(1.,1);
    auto t2 = std::make_shared<grid::atom::Type>(1.,2);
    auto r1 = std::make_shared<react>(t1,t2,t1,t2,1);
    auto r2 = std::make_shared<react>(t1,t2,t1,t2,1);

    g.AddReact(r1);
    g.AddReact(r2);
}

void react_standart() {
    Grid g(128.);

    auto t1 = make_shared<Type>(0,__COUNTER__,"O");
    auto t2 = make_shared<Type>(0,__COUNTER__,"X");
    auto t3 = make_shared<Type>(0,__COUNTER__,"Y");
    auto t4 = make_shared<Type>(0,__COUNTER__,"Z");

    grid::react::barrier p1(t1,t2,t3,t4,10,10,10);

    g.insert({1,2,3},make_shared<Atom>(t1));
    g.insert({1,3,3},make_shared<Atom>(t2));

    auto a1 = g.get({1,2,3});
    auto a2 = g.get({1,3,3});

    
    assert_eq(*a1->Material(),*t1);
    assert_eq(*a2->Material(),*t2);

    p1.Apply(a1,a2);

    grid::react::react r(t1,t2,t3,t4,10);


    assert_eq(*a1->Material(),*t3);
    assert_eq(*a2->Material(),*t4);
}

void grid_add_Atom() {
    Grid g(128.);

    auto t1 = std::make_shared<grid::atom::Type>(1.,1);
    auto t2 = std::make_shared<grid::atom::Type>(1.,2);

    auto a1 = std::make_shared<Atom>(t1);
    auto a2 = std::make_shared<Atom>(t2);
    auto a3 = std::make_shared<Atom>(t2);

    Vector p1 = Vector(50,50,50);
    Vector p2 = Vector(500,60,50);
    assert_eq(g.insert(p1,a1),grid::InsertionResults::OK);
    assert_eq(g.insert(p2,a2),grid::InsertionResults::OK);
    assert_eq(g.insert(p2,a3),grid::InsertionResults::RepeatedInsertion);

}

void grid_get_Atom() {
    Grid g(128.);

    auto t1 = std::make_shared<grid::atom::Type>(1.,1);
    auto t2 = std::make_shared<grid::atom::Type>(1.,2);

    auto a1 = std::make_shared<Atom>(t1);
    auto a2 = std::make_shared<Atom>(t2);
    auto a3 = std::make_shared<Atom>(t2);

    Vector p1 = Vector(50,50,50);
    Vector p2 = Vector(500,60,50);
    assert_eq(g.insert(p1,a1),grid::InsertionResults::OK);
    assert_eq(g.insert(p2,a2),grid::InsertionResults::OK);
    assert_eq(g.insert(p2,a3),grid::InsertionResults::RepeatedInsertion);

    auto aa1 = g.get(p1);
    auto aa2 = g.get(p2);

    assert_tst(a1==aa1);
    assert_tst(a2==aa2);

}

void grid_addLattice() {
    Grid g(10.);

    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);

    l.add(Vector(1,0,0),std::make_shared<Type>(1.,1.,"OMG"));

    g.AddLattice(Vector(0.,0,0),Vector(10.,10.,10),l);

    fmt::print("{}: atoms inserted: {}\n",__FUNCTION__,g.count());
}

void grid_cnt() {
    Grid g(128.);

    auto t1 = std::make_shared<grid::atom::Type>(1.,1);
    auto t2 = std::make_shared<grid::atom::Type>(1.,2);

    auto a1 = std::make_shared<Atom>(t1);
    auto a2 = std::make_shared<Atom>(t2);
    auto a3 = std::make_shared<Atom>(t2);

    Vector p1 = Vector(50,50,50);
    Vector p2 = Vector(500,60,50);

    g.insert(p1,a1);
    g.insert(p2,a2);

    assert_eq(g.count(),2);

}

void grid_add_lattive_ex() {
    Grid g(10.);


    auto t1 = std::make_shared<grid::atom::Type>(1.,1);

    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);

    g.AddLattice(Vector(0.,0,0),Vector(30.,30.,30),l);

    fmt::print("{}: atoms inserted: {}\n",__FUNCTION__,g.count());

}

void grid_clear_parallelep() {
    Grid g(10.);


    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);
    l.add({0,0,0},make_shared<Type>(0,0,"omg"));

    g.AddLattice(Vector(0.,0,0),Vector(10,10.,10),l);

    int cnt = g.count();


    g.ClearParallelep(Vector(5.,5.,5.),Vector(7.,7,7));


    assert_neq(cnt,g.count());



}

void grid_clear_parallelep_ex() {
    Grid g(10.);


    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);

    g.AddLattice(Vector(0.,0,0),Vector(30.,30.,30),l);


    g.ClearParallelep(Vector(0.,0.,0.),Vector(30.,30,30));

    assert_eq(g.count(),0);

}

void grid_erase() {
    Grid g(10.);



    auto t1 = std::make_shared<grid::atom::Type>(1.,1);
    auto a1 = std::make_shared<Atom>(t1);
    
    Vector p1 = Vector(100,100,100);
    Vector p2 = Vector(20,100,100);
    Vector p3 = Vector(20,20,20);


    assert_eq(g.count(),0);
    g.insert(p1,std::make_shared<Atom>(t1));
    assert_eq(g.count(),1);
    g.erase(p1);
    assert_eq(g.count(),0);
}

void grid_1m_insert_erase() {
    fmt::print("hello there!");
    Grid g(10.);



    auto t1 = std::make_shared<grid::atom::Type>(1.,1);
    auto a1 = std::make_shared<Atom>(t1);
    
    Vector p1 = Vector(100,100,100);
    Vector p2 = Vector(20,100,100);
    Vector p3 = Vector(20,20,20);

    Geometry geo(Vector(1.,0,0),Vector(0,1.,0),Vector(0,0,1));
    Lattice l(geo);

    l.add(Vector(1,1,1),std::make_shared<Type>(1,1));


    g.AddLattice(Vector(0.,0,0),Vector(10,10,10),l);

    int base_cnt = g.count();


    fmt::print("atoms: {}\n",base_cnt);


    assert_neq(base_cnt,0);


    for(size_t a = 0; a<1000;a++) {
        g.insert(p1,std::make_shared<Atom>(t1));
        assert_eq(g.count(),base_cnt+1);
        g.insert(p2,std::make_shared<Atom>(t1));
        assert_eq(g.count(),base_cnt+2);
        g.insert(p3,std::make_shared<Atom>(t1));
        assert_eq(g.count(),base_cnt+3);
        g.erase(p1);
        g.erase(p2);
        g.erase(p3);
    }
    assert_eq(g.count(),base_cnt);

    fmt::print("atoms: {}\n",base_cnt);
}

void grid_save() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"O")));


    fmt::print("grid: \n{}\n",g.to_xyz());

}

void grid_save_load() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"O")));



    auto gs = g.to_xyz();

    fmt::print("grid: \n{}\n",g.to_xyz());

    Grid g1(10.);

    g1.AddType(t);

    g1.from_xyz(gs);

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

    Grid g1(10.);

    g1.AddType(t1);
    g1.AddType(t2);

    g1.from_xyz(string_g);

    assert_eq(g.count(),g1.count());
}

void grid_radius_iterator() {
    Grid g(10.);
    fmt::print("1\n");
    auto t = std::make_shared<Type>(1,1,"O");
    fmt::print("1\n");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(t));

    fmt::print("1\n");

    Vector v(0,0,0);
    size_t cnt = 0;
    fmt::print("1\n");

    g.for_each(v,100,[t,&cnt](const Vector&,std::shared_ptr<Atom> atom) mutable
    {
        assert_tst((atom->Material())==t);
        cnt++;
    });
    assert_eq(cnt,1);
    fmt::print("3\n");


}

void grid_radius_iterator_exchecks() {
    Grid g(10.);
    fmt::print("1\n");
    auto t = std::make_shared<Type>(1,1,"O");
    fmt::print("1\n");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(10.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(100.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(101.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(101.,5.,1.),std::make_shared<Atom>(t));

    fmt::print("1\n");

    Vector v(100,0,0);
    size_t cnt = 0;
    fmt::print("1\n");


    g.for_each(v,10,[t,&cnt](const Vector&,std::shared_ptr<Atom> atom) mutable
    {
        assert_tst((atom->Material())==t);
        cnt++;
    });
    assert_eq(cnt,3);
    fmt::print("3\n");


}


void grid_radius_iterator_horizontal() {
    Grid g(10.);
    fmt::print("1\n");
    auto t = std::make_shared<Type>(1,1,"O");
    fmt::print("1\n");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(10.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(0.,100.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(0.,100.,2.),std::make_shared<Atom>(t));
    g.insert(Vector(0.,5.,1.),std::make_shared<Atom>(t));

    fmt::print("1\n");

    Vector v(0,100,0);
    size_t cnt = 0;
    fmt::print("1\n");

    g.for_each(v,10,[t,&cnt](const Vector&,std::shared_ptr<Atom> atom) mutable
    {
        assert_tst((atom->Material())==t);
        cnt++;
    });
    assert_eq(cnt,2);
    fmt::print("3\n");


}

void grid_radius_iterator_horizontal_z() {
    Grid g(10.);
    fmt::print("1\n");
    auto t = std::make_shared<Type>(1,1,"O");
    fmt::print("1\n");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(10.,1.,1.),std::make_shared<Atom>(t));
    g.insert(Vector(0.,0.,100.),std::make_shared<Atom>(t));
    g.insert(Vector(0.,1.,100.),std::make_shared<Atom>(t));
    g.insert(Vector(0.,5.,1.),std::make_shared<Atom>(t));

    fmt::print("1\n");

    Vector v(0,0,100);
    size_t cnt = 0;
    fmt::print("1\n");

    g.for_each(v,100,[t,&cnt](const Vector&,std::shared_ptr<Atom> atom) mutable
    {
        assert_tst((atom->Material())==t);
        cnt++;
    });
    assert_eq(cnt,5);
    fmt::print("3\n");


}

void grid_radius_iterator_simple() {
    Grid g(10.);

    auto t=make_shared<Type>(1.,1.,"O");

    g.insert(Vector(5.,5.,5.),std::make_shared<Atom>(t));


    Vector v(5,5,5);
    size_t cnt = 0;
    g.for_each(v,3,[t,&cnt](const Vector&,std::shared_ptr<Atom> atom) mutable
    {
        assert_tst(atom->Material()==t);
        cnt+=1;
    });
    assert_eq(cnt,1);
}

void grid_radius_iter_filtered() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    g.insert(Vector(1.,1.,1.),std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"O")));
    g.insert(Vector(1000.,1.,1.),std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"O")));


    Vector v(0,0,0);
    size_t cnt = 0;
    g.for_each(v,100,[&cnt](const Vector& pos,std::shared_ptr<Atom> atom) mutable
    {
        assert_tst(pos==Vector(1,1,1));
        cnt++;
    });
    assert_eq(cnt,1);


}

void grid_radius_iterator_ex() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    size_t b_cnt = 0;
    for (double d = 1; d<100; d+=1) {
        if(g.insert(Vector(d,1.,1.),std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"O")))==grid::InsertionResults::OK) b_cnt++;
    }


    Vector v(0,0,0);
    size_t cnt = 0;
    g.for_each(v,150,[&cnt](const Vector& pos,std::shared_ptr<Atom> atom) mutable
    {
        cnt++;
    });
    assert_eq(cnt,b_cnt);


}

void grid_radius_iterator_cyclic() {
    Grid g(10.);

    auto t = std::make_shared<Type>(1,1,"O");

    if(g.insert(Vector(1,1.,0.05),std::make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"O")))==grid::InsertionResults::OK);

    g.setPeriod(Vector(1.,1.,1));
    g.Cyclic<'x'>(true);
    g.Cyclic<'y'>(true);

    

    Vector v(0,0,0);
    size_t cnt = 0;
    g.for_each(v,0.06,[&cnt](const Vector& pos,std::shared_ptr<Atom> atom) mutable
    {
        cnt++;
        std::cout<<atom<<std::endl;
    });
    assert_eq(cnt,1)


}

void grid_cyclic_set() {
    Grid g(10.);

    g.Cyclic<'x'>(true);
    g.Cyclic<'y'>(true);
    g.Cyclic<'z'>(false);
}

void grid_cyclic_set_get() {
    Grid g(10.);

    g.Cyclic<'x'>(true);
    g.Cyclic<'y'>(true);
    g.Cyclic<'z'>(false);

    assert_eq(g.Cyclic<'x'>(),true);
    assert_eq(g.Cyclic<'y'>(),true);
    assert_eq(g.Cyclic<'z'>(),false);
}
using std::make_shared;

void grid_field_equal() {
    Grid g(10.);


    g.insert(Vector(0,0,0),make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"OOO")));
    g.insert(Vector(0,0,15),make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"OOO")));

    field::Equal Eq(1,0,10);

    Eq.Apply(g);

    assert_eq(g.get({0,0,0})->U(),1);
    assert_eq(g.get({0,0,15})->U(),0);
}


void grid_field_condenser() {
    Grid g(10.);


    g.insert({0,0,0},make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"OOO")));
    g.insert({0,0,10},make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"OOO")));
    g.insert({0,0,15},make_shared<Atom>(std::make_shared<grid::atom::Type>(1.,1.,"OOO")));

    field::ZCondenser cond(1,0,10);

    cond.Apply(g);

    assert_eq(g.get({0,0,0})->U(),0);
    assert_eq(g.get({0,0,10})->U(),1);
    assert_eq(g.get({0,0,15})->U(),0);
}

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

   g.for_each([](const Vector& pos,std::shared_ptr<Atom> atom) mutable{});

    exit(0);


}

void grid_ewald_hack_simple() {
    
    Grid g(10.);

    auto t1 = make_shared<Type>(1.,1,"O");

    g.insert({0,0,0},make_shared<Atom>(t1));

    field::ewald_hack hh(g);



    field::ZCondenser cond(1,0,10);

    cond.Apply(g);

    hh.calc_all();

}

void grid_check_reacts_calc() {
    
    Grid g(10.);

    auto t1 = make_shared<Type>(1.,1,"O");

    g.insert({0,0,0},make_shared<Atom>(t1));

    field::ewald_hack hh(g);



    field::ZCondenser cond(1,0,10);

    cond.Apply(g);

    hh.calc_all();

}

void grid_check_min_dist() {
    
    Grid g(10.);

    g.setPeriod({100,100,100});
    g.Cyclic<'x'>(true);
    g.Cyclic<'y'>(true);
    g.Cyclic<'z'>(true);

    auto t1 = make_shared<Type>(1.,1,"O");

    assert_eq(g.getMinDist({99,99,99},{1,1,1}),Vector(2,2,2));
    assert_eq(g.getMinDist({1,1,1},{99,99,99}),Vector(-2,-2,-2));
    assert_eq(g.getMinDist({1,1,1},{2,2,2}),Vector(1,1,1));
    assert_eq(g.getMinDist({2,2,2},{1,1,1}),Vector(-1,-1,-1));
    assert_eq(g.getMinDist({99,0,0},{1,0,0}),Vector(2,0,0));
    assert_eq(g.getMinDist({1,0,0},{9,0,0}),Vector(-2,0,0));

}

void gen_result() {

    Grid g(10*sgs::ANGSTROM);
    auto TElectrode = std::make_shared<Type>(0, __COUNTER__, "El",0);
    auto TElectrodeR = std::make_shared<Type>(0, __COUNTER__, "Elr",0);
    auto TElectrodeL = std::make_shared<Type>(0, __COUNTER__, "Ell",0);
    auto Oxygen = std::make_shared<Type>(-1* sgs::ELCHARGE, __COUNTER__, "O",6.4*powf(sgs::ANGSTROM,3));
    auto Oxygen_Intersittal = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "OI",6.4*powf(sgs::ANGSTROM,3));
    auto Hafnium = std::make_shared<Type>(2* sgs::ELCHARGE, __COUNTER__, "Hf",21.88*powf(sgs::ANGSTROM,3));
    auto OxygenVacancy_Neutral = std::make_shared<Type>(0, __COUNTER__, "Vo",6.4*powf(sgs::ANGSTROM,3));
    auto IntersitialPosition = std::make_shared<Type>(0, __COUNTER__, "Ip",6.4*powf(sgs::ANGSTROM,3));
    auto OxygenVacancy_Charged = std::make_shared<Type>(-1 * sgs::ELCHARGE, __COUNTER__, "Vo",6.4*powf(sgs::ANGSTROM,3));

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
    lHfO2.add({0.932151, 0.330122, 0.652906}, OxygenVacancy_Neutral);
    lHfO2.add({0.067849, 0.830122, 0.847094}, OxygenVacancy_Neutral);
    lHfO2.add({0.067849, 0.669878, 0.347094}, OxygenVacancy_Neutral);
    lHfO2.add({0.932151, 0.169878, 0.152906}, OxygenVacancy_Neutral);
    lHfO2.add({0.78284, 0.43966, 0.7821}, Oxygen_Intersittal);
    lHfO2.add({0.21066, 0.5392, 0.22844}, Oxygen_Intersittal);
    lHfO2.add({0.61202, 0.9786, 0.39605}, Oxygen_Intersittal);
    lHfO2.add({0.38798, 0.0213, 0.60394}, Oxygen_Intersittal);

    double size_x = 30 * sgs::ANGSTROM;
    double size_y = 30 * sgs::ANGSTROM;
    double size_z = 30 * sgs::ANGSTROM;

    double dist_electrode = 0.3 * sgs::ANGSTROM;

    double electrode_end = 1 * sgs::ANGSTROM;

    double oxyde_begin = electrode_end + dist_electrode;

    double oxyde_end = 20* sgs::ANGSTROM + oxyde_begin;


    double el_begin = oxyde_end+dist_electrode;


    double electrode2_end = 3 * sgs::ANGSTROM + el_begin;

    g.AddLattice({0., 0, 0}, Vector(size_x, size_y, electrode_end), lElectrode);
    g.AddLattice({0, 0, oxyde_begin}, Vector(size_x, size_y, oxyde_end), lHfO2);
    g.AddLattice({0, 0, el_begin}, Vector(size_x, size_y, electrode2_end), lElectrode);

    double size_electrode = sgs::ANGSTROM * 10;
    double size_electrode_z = sgs::ANGSTROM * 10;
    double electrode_begin_xy = sgs::ANGSTROM * 10;
    double electrode_begin_z = electrode_end;

    Vector ElCube_s(electrode_begin_xy, electrode_begin_xy, electrode_begin_z);
    Vector ElCube_e(electrode_begin_xy + size_electrode, electrode_begin_xy + size_electrode, electrode_begin_z + size_electrode_z);

    g.ClearParallelep(ElCube_s, ElCube_e);

    int cnt = g.count();

    g.AddLattice(ElCube_s, ElCube_e, lElectrode);

    double struct_end = electrode2_end;

    {
        std::ofstream out("generated_result.xyz");

        out << g.to_xyz(1 / sgs::ANGSTROM);

        out.close();
    }
}

#define add_test(fn) tests[#fn] = fn

std::unordered_map<std::string,std::function<void()>> tests;

void init_tests() {
    add_test(all_work);
    add_test(math_modulo_test_double);
    add_test(test_basic_sgs);
    add_test(test_vector_creation);
    add_test(vector_add);
    add_test(vector_cmp_check);
    add_test(vector_hash);
    add_test(vector_cmp);
    add_test(test_must_crush);
    add_test(geometry_create);
    add_test(geometry_to);
    add_test(geometry_from);
    add_test(cubic_chunk_creation);
    add_test(cubic_chunk_insert);
    add_test(cubic_chunk_multy);
    add_test(cubic_chunk_insert_get);
    add_test(cubic_chunk_iterate);
    add_test(cubic_chunk_count);
    add_test(cubic_chunk_erase);
    add_test(atom_type_creation_and_cmp);
    add_test(lattice_creation_check);
    add_test(lattice_creation_check_incorrect_add);
    add_test(check_atom_create);
    add_test(check_reaction_create);
    add_test(check_reaction_check_atoms);
    add_test(create_grid);
    add_test(grid_add_Reaction);
    add_test(grid_addLattice);
    add_test(grid_add_Atom);
    add_test(grid_get_Atom);
    add_test(grid_cnt);
    add_test(grid_add_lattive_ex);
    add_test(grid_clear_parallelep);
    add_test(grid_clear_parallelep_ex);
    add_test(grid_erase);
    add_test(grid_1m_insert_erase);
    add_test(grid_save);
    add_test(grid_save_load);
    add_test(grid_save_load_big);
    add_test(grid_radius_iterator);
    add_test(grid_radius_iterator_simple);
    add_test(grid_radius_iter_filtered);
    add_test(grid_radius_iterator_ex);
    add_test(grid_radius_iterator_exchecks);
    add_test(grid_radius_iterator_horizontal);
    add_test(grid_radius_iterator_horizontal_z);
    add_test(grid_cyclic_set);
    add_test(grid_cyclic_set_get);
    add_test(grid_radius_iterator_cyclic);
    add_test(grid_field_equal);
    add_test(grid_field_condenser);
    add_test(grid_ewald_hack_simple_demo);
    add_test(react_standart);
    add_test(grid_check_min_dist);
    add_test(gen_result);
}

bool run_test(string s) {
    auto itest = tests.find(s);
    if(tests.end()!=itest) {
        auto& [name,callback] = *itest;
        
        fmt::print("running test '{}'\n",s);
        
        try {
            callback();
        }
        catch (std::exception e) {
            fmt::print("exception happened: '{}'\n",e.what());
            return false;
        }
        fmt::print("runned successfully\n");
        return true;
    }
    else
    {
        fmt::print("test: '{}' not found\n",s);
        return false;
    }
    return false;
}



int main(int argc, char** argv) {
    if(argc<=1) { std::printf("incorrect amount of arguments: {}",argc); return -1;}
    auto test = string(argv[1]);
    size_t t_cnt = 0;
    
    init_tests();

    if(run_test(test)) {
        return 0;
    }
    return -1;
}