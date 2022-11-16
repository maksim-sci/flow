#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <memory>
#include <unordered_map>

#include <geometry/Vector.hpp>

#include <math/modulo.hpp>

#include "cubicchunk.hpp"
#include "atom/atom.hpp"
#include "atom/type.hpp"
#include "react/react.hpp"
#include "lattice.hpp"

namespace grid
{

    using std::string;
    struct LatticePos
    {
        geometry::Vector llim{};
        geometry::Vector rlim{};

        grid::Lattice lat{};
    };

    typedef std::unordered_map<Vector, std::shared_ptr<atom::Atom>> AtomMap;
    typedef std::unordered_map<Vector, std::shared_ptr<chunk::CubicChunk>> ChunkMap;

    class Grid
    {
        ChunkMap chunks;

        geometry::Vector llim;
        geometry::Vector rlim;
        geometry::Vector delta;

        std::unordered_map<std::string, std::shared_ptr<atom::Type>> types;

        std::vector<std::shared_ptr<grid::react::React>> reacts;
        std::vector<LatticePos> lattices;

        double size_chunk;

        bool cycle_x, cycle_y, cycle_z;

    public:

        geometry::Vector calcChunkPos(const geometry::Vector &pos) const;

        class GridIteratorOnceSinglePass
        {

            grid::ChunkMap::const_iterator citer;
            grid::ChunkMap::const_iterator citerend;
            grid::AtomMap::const_iterator aiterend;

        public:
            grid::AtomMap::const_iterator aiter;

            GridIteratorOnceSinglePass(const Grid &g);
            inline bool operator==(const GridIteratorOnceSinglePass &g) const;
            GridIteratorOnceSinglePass &operator++();


            inline bool operator!=(const Grid::GridIteratorOnceSinglePass &g) const {return g!=*this;};


            void end();
        };

        class GridIteratorDistLimSinglePass
        {
            double distance;
            double chunk_size;

            Vector pos;
            Vector c0;

            Grid *grid;

            int translation_x, translation_y, translation_z;
            int mtranslation_x, mtranslation_y, mtranslation_z;

            Vector translations;

            ChunkMap *chunks;

            size_t chunksAround;

            grid::AtomMap::const_iterator aiterend;

            bool finished;

            double x, y, z;
            double mx, my, mz;

            void iterateTillCorrectOrEnd();

        public:
            grid::AtomMap::const_iterator aiter;

            GridIteratorDistLimSinglePass(Grid &g, double dist, const Vector &p);
            bool operator==(const GridIteratorDistLimSinglePass &g) const;
            GridIteratorDistLimSinglePass &operator++();

            inline bool Finished() { return finished; };
        };

    public:
        inline Grid(double chunkSize) : llim(0, 0, 0), rlim(0, 0, 0), reacts(), lattices(), size_chunk(chunkSize), cycle_x(false), cycle_y(false), cycle_z(false),delta(rlim-llim){};

        inline void AddReact(std::shared_ptr<grid::react::React> r) { reacts.push_back(r); };
        inline void AddType(std::shared_ptr<Type> t) { types[t->Name()] = t; };
        void AddLattice(const geometry::Vector &l, const geometry::Vector &r, grid::Lattice latt);
        void ClearParallelep(const geometry::Vector &l, const geometry::Vector &r);

        inline const auto &Reacts() const;
        inline const auto &Lattices() const;
        inline const auto &Llim() const;
        inline const auto &Rlim() const;
        double Chunk_Size() const;
        size_t count() const;

        inline geometry::Vector Sizes() const { return rlim - llim; };

        InsertionResults insert(const geometry::Vector &pos, std::shared_ptr<atom::Atom> a);
        const std::shared_ptr<atom::Atom> get(const geometry::Vector &pos) const;
        bool erase(const geometry::Vector &p);

        inline const ChunkMap &Chunks() const { return chunks; };

        inline auto begin() const { return GridIteratorOnceSinglePass(*this); };
        inline auto end() const
        {
            auto a = GridIteratorOnceSinglePass(*this);
            a.end();
            return a;
        }
        inline auto beginFilterDistance(double dist, const Vector &p) { return GridIteratorDistLimSinglePass(*this, dist, p); }; // WARNING!!!! UNSAFE
        inline void setPeriod(const Vector &v) { rlim = llim + v; };

        template <char coord>
        inline bool Cyclic() const
        {
            static_assert(coord == 'x' || coord == 'y' || coord == 'z');
            if constexpr (coord == 'x')
                return cycle_x;
            if constexpr (coord == 'y')
                return cycle_y;
            if constexpr (coord == 'z')
                return cycle_z;
        }

        template <char coord>
        inline void Cyclic(bool cyclic)
        {
            static_assert(coord == 'x' || coord == 'y' || coord == 'z');
            if constexpr (coord == 'x')
                cycle_x = cyclic;
            if constexpr (coord == 'y')
                cycle_y = cyclic;
            if constexpr (coord == 'z')
                cycle_z = cyclic;
        }

        Vector getCycledVector(const Vector &v) const;

        string to_xyz() const;
        string to_xyz(double mult) const;
        void from_xyz(const string &s);
        void from_xyz(const string &s, double div);


        inline Vector getMinDist(const Vector& a, const Vector& b) {
            auto p1 = getCycledVector(b);
            Vector dist = (p1-a-llim);
            auto& [dx,dy,dz] = dist;
            auto& [sx,sy,sz] = delta;
            if(cycle_x) {
                dx = std::abs(dx);
                if(dx>sx/2) {
                    dx = sx-dx;
                }
            }
            if(cycle_y) {
                dy = std::abs(dy);
                if(dy>sy/2) {
                    dy = sy-dy;
                }
            }
            if(cycle_z) {
                dz = std::abs(dz);
                if(dx>sz/2) {
                    dz = sz-dz;
                }
            }

            return Vector(dx,dy,dz);
        };
    };

    inline const auto &grid::Grid::Reacts() const { return reacts; };

    inline const auto &Grid::Lattices() const { return lattices; };

    inline const auto &Grid::Llim() const { return llim; };

    inline const auto &Grid::Rlim() const { return rlim; };

    inline bool Grid::GridIteratorOnceSinglePass::operator==(const Grid::GridIteratorOnceSinglePass &g) const
    
    {
    
        if (aiter == aiterend && citer == citerend)
    
            return true;
    
        return false;
    
    };

}