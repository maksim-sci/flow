


#include <iostream>
#include <memory>
#include <string>
#include <map>
#include <fmt/format.h>
#include <algorithm>
#include <unordered_map>
#include <vector>

#include <sgs.hpp>
#include "fmt/core.h"
#include "macro_switch.hpp"
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
#include <grid/grid.hpp>
#include <math/modulo.hpp>
#include <field/condenser.hpp>
#include <field/equal.hpp>
#include <field/ewald_hack.hpp>
#include <algo/kmk.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>

#include <grid/react/barrier.hpp>

#include <assertions.h>

#include <chrono>

template <class Tp>
[[gnu::always_inline]] inline void DoNotOptimize(Tp const& value) {
  asm volatile("" : : "r,m"(value) : "memory");
}

template <class Tp>
[[gnu::always_inline]] inline void DoNotOptimize(Tp& value) {
#if defined(__clang__)
  asm volatile("" : "+r,m"(value) : : "memory");
#else
  asm volatile("" : "+m,r"(value) : : "memory");
#endif
}

using std::string;
using geometry::Vector;
using geometry::Geometry;
using grid::chunk::CubicChunk;

using grid::atom::Atom;
using grid::atom::Type;

using std::shared_ptr;

using std::make_shared;

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

#define add_test(fn) tests[#fn] = fn

struct pos_atom {
    Vector pos;
    shared_ptr<Atom> atom;
};


void simple_check() {
    std::unordered_map<Vector, shared_ptr<Atom>> chunk_ht;
    std::vector<pos_atom> chunk_vector;


    auto t1=make_shared<Type>(0.,1,"OXYGEN");
    auto a1=make_shared<Atom>(t1);
    auto a2=make_shared<Atom>(t1);
    auto a3=make_shared<Atom>(t1);
    auto a4=make_shared<Atom>(t1);
    auto a5=make_shared<Atom>(t1);
    auto a6=make_shared<Atom>(t1);
    auto a7=make_shared<Atom>(t1);
    auto a8=make_shared<Atom>(t1);

    Vector p1{1, 1, 1};
    Vector p2{1, 1, 2};
    Vector p3{1, 2, 1};
    Vector p4{1, 3, 1}; 
    Vector p5{2, 1, 1};
    Vector p6{3, 1, 1};
    Vector p7{4, 1, 1};
    Vector p8{1, 2, 9};
    
    auto data = {pos_atom{p1,a1},pos_atom{p2,a2},pos_atom{p3,a3},pos_atom{p4,a4},pos_atom{p5,a5},pos_atom{p6,a6},pos_atom{p7,a7},pos_atom{p8,a8}};
    for(auto& [pos,atom]:data) {
        chunk_ht.insert({pos,atom});
        chunk_vector.push_back({pos,atom});
    }
}


void insert_check() {
    std::unordered_map<Vector, shared_ptr<Atom>> chunk_ht;
    std::vector<pos_atom> chunk_vector;


    
    constexpr size_t count = 100000;
    auto perf_check = [](auto& container,auto func,std::string msg) {
        auto t1=make_shared<Type>(0.,1,"OXYGEN");
        auto a1=make_shared<Atom>(t1);
        auto a2=make_shared<Atom>(t1);
        auto a3=make_shared<Atom>(t1);
        auto a4=make_shared<Atom>(t1);
        auto a5=make_shared<Atom>(t1);
        auto a6=make_shared<Atom>(t1);
        auto a7=make_shared<Atom>(t1);
        auto a8=make_shared<Atom>(t1);

        Vector p1{1, 1, 1};
        Vector p2{1, 1, 2};
        Vector p3{1, 2, 1};
        Vector p4{1, 3, 1}; 
        Vector p5{2, 1, 1};
        Vector p6{3, 1, 1};
        Vector p7{4, 1, 1};
        Vector p8{1, 2, 9};
        
        auto data = {pos_atom{p1,a1},pos_atom{p2,a2},pos_atom{p3,a3},pos_atom{p4,a4},pos_atom{p5,a5},pos_atom{p6,a6},pos_atom{p7,a7},pos_atom{p8,a8}};
        
        std::chrono::high_resolution_clock clock;
        auto t_begin = clock.now();
        for(size_t i = 1; i<count;i++) {
            container.clear();
            for(auto& [pos,atom]:data) {
                func(container,{pos,atom});
                DoNotOptimize(pos);
                DoNotOptimize(atom);
            }
        }
        auto t_end = clock.now();

        std::chrono::duration< double > fs = t_end  - t_begin;
        std::chrono::milliseconds d = std::chrono::duration_cast< std::chrono::milliseconds >( fs );
        std::cout<<msg<<" "<<d.count()<<std::endl;
    };

    constexpr auto vector_insert = [](auto& container,const pos_atom element){container.push_back(element);};
    constexpr auto ht_insert = [](auto& container,const std::pair<Vector,shared_ptr<Atom>> element){container.insert(element);};

    perf_check(chunk_ht,ht_insert,"ht: ");
    perf_check(chunk_vector,vector_insert,"vector: ");
    
}

void check_iteration() {
    std::unordered_map<Vector, shared_ptr<Atom>> chunk_ht;
    std::vector<pos_atom> chunk_vector;


    
    constexpr size_t count = 1000000;
    auto perf_check = [](auto& container,auto func,std::string msg) {
        auto t1=make_shared<Type>(0.,1,"OXYGEN");
        auto a1=make_shared<Atom>(t1);
        auto a2=make_shared<Atom>(t1);
        auto a3=make_shared<Atom>(t1);
        auto a4=make_shared<Atom>(t1);
        auto a5=make_shared<Atom>(t1);
        auto a6=make_shared<Atom>(t1);
        auto a7=make_shared<Atom>(t1);
        auto a8=make_shared<Atom>(t1);

        Vector p1{1, 1, 1};
        Vector p2{1, 1, 2};
        Vector p3{1, 2, 1};
        Vector p4{1, 3, 1}; 
        Vector p5{2, 1, 1};
        Vector p6{3, 1, 1};
        Vector p7{4, 1, 1};
        Vector p8{1, 2, 9};
        
        auto data = {pos_atom{p1,a1},pos_atom{p2,a2},pos_atom{p3,a3},pos_atom{p4,a4},pos_atom{p5,a5},pos_atom{p6,a6},pos_atom{p7,a7},pos_atom{p8,a8}};
        
        container.clear();
        for(auto& [pos,atom]:data) {
            func(container,{pos,atom});
        }
        std::chrono::high_resolution_clock clock;
        auto t_begin = clock.now();
        for(size_t i = 1; i<count;i++) {
            for(auto&[pos,atom]:container) {
                DoNotOptimize(pos);
                DoNotOptimize(atom);
            };
        }
        auto t_end = clock.now();

        std::chrono::duration< double > fs = t_end  - t_begin;
        std::chrono::milliseconds d = std::chrono::duration_cast< std::chrono::milliseconds >( fs );
        std::cout<<msg<<" "<<d.count()<<std::endl;
    };

    constexpr auto vector_insert = [](auto& container,const pos_atom element){container.push_back(element);};
    constexpr auto ht_insert = [](auto& container,const std::pair<Vector,shared_ptr<Atom>> element){container.insert(element);};

    perf_check(chunk_ht,ht_insert,"ht: ");
    perf_check(chunk_vector,vector_insert,"vector: ");
    
}

void check_lookup() {
    std::unordered_map<Vector, shared_ptr<Atom>> chunk_ht;
    std::vector<pos_atom> chunk_vector;


    
    constexpr size_t count = 1000000;
    auto perf_check = [](auto& container,auto inserter,auto checker,std::string msg) {
        auto t1=make_shared<Type>(0.,1,"OXYGEN");
        auto a1=make_shared<Atom>(t1);
        auto a2=make_shared<Atom>(t1);
        auto a3=make_shared<Atom>(t1);
        auto a4=make_shared<Atom>(t1);
        auto a5=make_shared<Atom>(t1);
        auto a6=make_shared<Atom>(t1);
        auto a7=make_shared<Atom>(t1);
        auto a8=make_shared<Atom>(t1);

        Vector p1{1, 1, 1};
        Vector p2{1, 1, 2};
        Vector p3{1, 2, 1};
        Vector p4{1, 3, 1}; 
        Vector p5{2, 1, 1};
        Vector p6{3, 1, 1};
        Vector p7{4, 1, 1};
        Vector p8{1, 2, 9};
        
        auto data = {pos_atom{p1,a1},pos_atom{p2,a2},pos_atom{p3,a3},pos_atom{p4,a4},pos_atom{p5,a5},pos_atom{p6,a6},pos_atom{p7,a7},pos_atom{p8,a8}};
        
        container.clear();
        for(auto& [pos,atom]:data) {
            inserter(container,{pos,atom});
        }
        std::chrono::high_resolution_clock clock;
        auto t_begin = clock.now();

        Vector pos1{9,9,9};
        Vector pos2{8,3,1};
        Vector pos3{1,5,1};

        
        for(size_t i = 1; i<count;i++) {
            DoNotOptimize(checker(container,pos1));
            DoNotOptimize(checker(container,p1));
            DoNotOptimize(checker(container,pos2));
            DoNotOptimize(checker(container,p2));
            DoNotOptimize(checker(container,pos3));
            DoNotOptimize(checker(container,p3));
        }

        auto t_end = clock.now();

        std::chrono::duration< double > fs = t_end  - t_begin;
        std::chrono::milliseconds d = std::chrono::duration_cast< std::chrono::milliseconds >( fs );
        std::cout<<msg<<" "<<d.count()<<std::endl;
    };

    constexpr auto vector_insert = [](auto& container,const pos_atom element){container.push_back(element);};
    constexpr auto ht_insert = [](auto& container,const std::pair<Vector,shared_ptr<Atom>> element){container.insert(element);};

    constexpr auto checker_ht = [](const auto& container, const auto& pos)->shared_ptr<Atom>{
        auto iter = container.find(pos);
        if(iter!=container.end()) {
            return iter->second;
        }
        return nullptr;
    };

    constexpr auto checker_vector = [](const auto& container, const auto& pos)->shared_ptr<Atom>{
        for(const auto& [pos1,atom]:container) {
            if(pos==pos1) return atom;
        }
        return nullptr;
    };

    perf_check(chunk_ht,ht_insert,checker_ht,"ht: ");
    perf_check(chunk_vector,vector_insert,checker_vector,"vector: ");
    
}

std::unordered_map<std::string,std::function<void()>> tests;

void init_tests() {
    add_test(insert_check);
    add_test(check_iteration);
    add_test(check_lookup);
    
    
    

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
    if(argc<=1) { std::printf("incorrect amount of arguments: %d",argc); return -1;}
    auto test = string(argv[1]);
    size_t t_cnt = 0;
    
    init_tests();

    if(run_test(test)) {
        return 0;
    }
    return -1;
}