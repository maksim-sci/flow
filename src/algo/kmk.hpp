#pragma once

#include "grid/lattice.hpp"
#include <list>
#include <memory>
#include <unordered_map>
#include <random>
#include <vector>

#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/react/react.hpp>
#include <grid/grid.hpp>

namespace algo {

class kmk {
    public:
    using ptr_react = std::shared_ptr<grid::react::react>;
    bool cache{false};

    struct react_data {
        std::shared_ptr<grid::atom::Atom> f;
        std::shared_ptr<grid::atom::Atom> s;
        geometry::Vector fp;
        geometry::Vector sp;
        ptr_react r;
        double chance;
    };

    struct cache_elem {
        std::shared_ptr<grid::atom::Atom> f;
        std::shared_ptr<grid::atom::Atom> s;
        ptr_react r;
        geometry::Vector fp;
        geometry::Vector sp;
    };

    std::vector<cache_elem> cacheData{};

    grid::Grid* g;
    double sum{0};
    double time{0};
    bool apply{false};

    std::vector<ptr_react> reacts;

    std::unordered_multimap<std::shared_ptr<grid::atom::Atom>, react_data> calculated_reacts;
    std::unordered_multimap<std::shared_ptr<grid::atom::Atom>, react_data> recieved_reacts;

    void add(ptr_react r);

    void remove(ptr_react r);

    void recalc();

    std::pair<bool, react_data> chooseReact();

    void Cache(bool mode);

    void processReact(react_data &data);

    size_t Count() const;
    double Sum() const;
    double Time() const;
    void ClearTime();

    std::pair<bool, react_data> findAndProcessReact();

    kmk(grid::Grid *_g);
};

}