


#include <iostream>
#include <string>
#include <map>
#include <fmt/format.h>

#include <sgs.hpp>
#include "macro_switch.hpp"
#include <geometry/vector.hpp>
#include <geometry/geometry.hpp>
#include <grid/cubicchunk.hpp>
#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/lattice.hpp>
#include <grid/react/react.hpp>
#include <grid/react/standart.hpp>
#include <grid/grid.hpp>
#include <math/modulo.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>
#include <field/ewald_hack.hpp>

#include <grid/react/standart.hpp>

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


    assert_eq(a.insert({0.,0.,0.},std::make_shared<Atom>(t)),grid::InsertionResults::OK);
    assert_eq(a.insert({.5,0.,0.},std::make_shared<Atom>(t)),grid::InsertionResults::OK);
    assert_eq(a.insert({0,.5,0.},std::make_shared<Atom>(t)),grid::InsertionResults::OK);
    assert_eq(a.insert({0,.5,0.},std::make_shared<Atom>(t)),grid::InsertionResults::RepeatedInsertion);

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

void react_puasson() {
    Grid g(128.);

    auto t1 = make_shared<Type>(0,__COUNTER__,"O");
    auto t2 = make_shared<Type>(0,__COUNTER__,"X");
    auto t3 = make_shared<Type>(0,__COUNTER__,"Y");
    auto t4 = make_shared<Type>(0,__COUNTER__,"Z");

    grid::react::standart p1(t1,t2,t3,t4,10,10,10);

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


#define Test_case(fn,code) CASE(code): {fn();t_cnt++;break;}

int main(int argc, char** argv) {
    if(argc<=1) { std::printf("incorrect amount of arguments: {}",argc); return -1;}
    auto arg1 = string(argv[1]);
    std::cout<<"running test: "<<arg1<<std::endl;
    size_t t_cnt = 0;
    SWITCH (arg1) {
    Test_case(all_work,"aw");
    Test_case(math_modulo_test_double,"mmtd");
    Test_case(test_basic_sgs,"bsgs");
    Test_case(test_vector_creation,"vc");
    Test_case(vector_add,"va");
    Test_case(vector_cmp_check,"vcc");
    Test_case(vector_hash,"vhash");
    Test_case(vector_cmp,"vcmp");
    Test_case(test_must_crush,"mf");
    Test_case(geometry_create,"gt");
    Test_case(geometry_to,"gtto");
    Test_case(geometry_from,"gtfrom");
    Test_case(cubic_chunk_creation,"ccc");
    Test_case(cubic_chunk_insert,"cci");
    Test_case(cubic_chunk_multy,"ccm");
    Test_case(cubic_chunk_insert_get,"ccig");
    Test_case(cubic_chunk_iterate,"cciter");
    Test_case(cubic_chunk_count,"ccc1");
    Test_case(cubic_chunk_erase,"cce");
    Test_case(atom_type_creation_and_cmp,"atcr");
    Test_case(lattice_creation_check,"lcc");
    Test_case(lattice_creation_check_incorrect_add,"lccia");
    Test_case(check_atom_create,"cac");
    Test_case(check_reaction_create,"crc");
    Test_case(check_reaction_check_atoms,"crca");
    Test_case(create_grid,"cg");
    Test_case(grid_add_Reaction,"gar");
    Test_case(grid_addLattice,"gal");
    Test_case(grid_add_Atom,"gaa");
    Test_case(grid_get_Atom,"gga");
    Test_case(grid_cnt,"gc");
    Test_case(grid_add_lattive_ex,"gale");
    Test_case(grid_clear_parallelep,"gcp");
    Test_case(grid_clear_parallelep_ex,"gcpe");
    Test_case(grid_erase,"ge");
    Test_case(grid_1m_insert_erase,"g1mie");
    Test_case(grid_save,"gs");
    Test_case(grid_save_load,"gsl");
    Test_case(grid_save_load_big,"gslb");
    Test_case(grid_radius_iterator,"gri");
    Test_case(grid_radius_iterator_simple,"gris");
    Test_case(grid_radius_iter_filtered,"grif");
    Test_case(grid_radius_iterator_ex,"grie");
    Test_case(grid_radius_iterator_exchecks,"griecss");
    Test_case(grid_radius_iterator_horizontal,"grieh");
    Test_case(grid_radius_iterator_horizontal_z,"griehz");
    Test_case(grid_cyclic_set,"gcs");
    Test_case(grid_cyclic_set_get,"gcsg");
    Test_case(grid_radius_iterator_cyclic,"gric");
    Test_case(grid_field_equal,"gfe");
    Test_case(grid_field_condenser,"gfc");
    Test_case(grid_ewald_hack_simple_demo,"gehsd");
    Test_case(react_puasson,"rp");
    }

    if(t_cnt==1) {printf("done\n");return 0;}
    return -1;
}