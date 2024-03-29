#pragma once

#include <functional>

#include "../geometry/vector.hpp"
#include <unordered_map>

#include <memory>
#include "atom/atom.hpp"


namespace grid
{
    using for_each_callbak = std::function<void(const geometry::Vector&,std::shared_ptr<atom::Atom>&)>;
    using for_each_const_callbak = std::function<void(const geometry::Vector&,std::shared_ptr<atom::Atom>&)>;//TODO shared ptr must be const

    enum InsertionResults
    {
        OK,
        RepeatedInsertion,
        InsertionOutOfBounds,
        UNKNOWN
    };

    namespace chunk
    {

        using geometry::Vector;
        using std::unordered_map;
        class CubicChunk
        {

            Vector pos;
            double size;

            unordered_map<Vector, std::shared_ptr<atom::Atom>> atoms;

        public:
            CubicChunk();
            CubicChunk(double x, double y, double z, double size);
            CubicChunk(const Vector &pos, double size);

            size_t cnt() const;

            Vector getLeftLimit() const;
            Vector getRightLimit() const;
            double dim_size() const;
            bool inPosIn(const Vector &v) const;
            inline auto begin() const;
            inline auto end() const;
            size_t count() const;

            std::shared_ptr<atom::Atom> get(const Vector &pos) const;

            InsertionResults insert(const Vector& pos, std::shared_ptr<atom::Atom>& atom);
            bool erase(const Vector &pos);

            void for_each(for_each_const_callbak callback) const;
        };

        inline auto CubicChunk::begin() const { return atoms.begin(); };

        inline auto CubicChunk::end() const { return atoms.end(); };
    }

    

}