#pragma once

#include "grid/lattice.hpp"
#include <list>
#include <memory>
#include <unordered_map>
#include <random>

#include <grid/atom/type.hpp>
#include <grid/atom/atom.hpp>
#include <grid/react/react.hpp>
#include <grid/grid.hpp>

namespace algo {

class kmk {
    public:
    using ptr_react = std::shared_ptr<grid::react::react>;

    struct react_data {
        std::shared_ptr<grid::atom::Atom> f;
        std::shared_ptr<grid::atom::Atom> s;
        geometry::Vector fp;
        geometry::Vector sp;
        ptr_react r;
        double chance;
    };

    grid::Grid* g;
    double sum{0};
    double time{0};

    std::vector<ptr_react> reacts;

    std::unordered_multimap<std::shared_ptr<grid::atom::Atom>, react_data> calculated_reacts;
    std::unordered_multimap<std::shared_ptr<grid::atom::Atom>, react_data> recieved_reacts;

    void add(ptr_react r);

    void remove(ptr_react r);

    void recalc();

    std::pair<bool, react_data> chooseReact();

    void processReact(react_data &data);

    size_t Count() const;
    double Sum() const;
    double Time() const;
    void ClearTime();

    std::pair<bool, react_data> findAndProcessReact();

    kmk(grid::Grid *_g);
};

}