#pragma once

#include <map>
#include <set>

#include "../geometry/vector.hpp"
#include "../geometry/geometry.hpp"
#include "atom/type.hpp"

namespace grid
{
    using geometry::Geometry;
    using geometry::Vector;
    using grid::atom::Type;
    using grid::atom::TypeId;
    using std::map;
    using std::set;

    class Lattice
    {

        set<std::shared_ptr<Type>> types;

        std::map<Vector, std::shared_ptr<Type>> positions;

        Geometry translations;

    public:
        inline Lattice() : translations(), types(), positions(){};
        inline Lattice(Geometry _translations) : translations(_translations), types(), positions(){};

        void add(Vector pos, std::shared_ptr<Type> atype);

        inline bool hasType(std::shared_ptr<Type> type) const { return (types.find(type) != types.end()); };

        inline auto begin() { return positions.begin(); };
        inline auto end() { return positions.end(); };

        inline const Geometry &Translations() const { return translations; };

        inline const set<std::shared_ptr<Type>> &Types() const { return types; };
        inline const std::map<Vector, std::shared_ptr<Type>> &Positions() const { return positions; };
    };
}