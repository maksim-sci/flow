#include <cmath>
#include <memory>
#include <sstream>

#include "grid.hpp"
#include "atom/atom.hpp"
#include "react/react.hpp"

namespace grid
{
    using geometry::Vector;

    Vector Grid::calcChunkPos(const Vector &pos) const // has limited accuracy!
    {
        auto vpos = (pos - llim);

        vpos /= size_chunk;

        constexpr double assimp_const(1 + 1e-10);
        double x = std::floor(vpos.x * assimp_const) * size_chunk;
        double y = std::floor(vpos.y * assimp_const) * size_chunk;
        double z = std::floor(vpos.z * assimp_const) * size_chunk;
        return geometry::Vector(x, y, z);
    }

    InsertionResults Grid::insert(const Vector &pos, std::shared_ptr<atom::Atom> a)
    {
        auto vpos = this->calcChunkPos(pos);

        auto pchunk = chunks.find(vpos);
        if (pchunk == chunks.end())
        {
            chunks[vpos] = std::make_shared<grid::chunk::CubicChunk>(vpos, size_chunk);
            pchunk = chunks.find(vpos);
            if (pchunk == chunks.end())
            {
                return InsertionResults::UNKNOWN;
            }
        }

        return pchunk->second->insert(pos, a);
    };
    const std::shared_ptr<atom::Atom> Grid::get(const Vector &pos) const
    {
        auto vpos = this->calcChunkPos(pos);

        auto pchunk = chunks.find(vpos);

        if (pchunk == chunks.end())
            return nullptr;

        return pchunk->second->get(pos);
    };

    bool Grid::erase(const geometry::Vector &p)
    
    {
    
        auto pc = calcChunkPos(p);
    
        auto chunk = chunks.find(p);
    
        if (chunk == chunks.end())
    
            return false;
    
        return chunk->second->erase(p);
    
    };

    void Grid::AddLattice(const geometry::Vector &l, const geometry::Vector &r, grid::Lattice latt)
    {
        lattices.push_back(LatticePos{l, r, latt});

        Vector latlim = llim + r;

        if (latlim.x > rlim.x)
        {
            rlim.x = latlim.x;
        };
        if (latlim.y > rlim.y)
        {
            rlim.y = latlim.y;
        };
        if (latlim.z > rlim.z)
        {
            rlim.z = latlim.z;
        };

        for (auto &a : latt.Types())
        {
            AddType(a);
        }
        Vector p0 = l + llim;
        Vector p1 = r + llim;

        const auto translations = latt.Translations();

        const auto translation1 = translations.A();
        const auto translation2 = translations.B();
        const auto translation3 = translations.C();

        // TODO fix this!!!!
        // works only with good translations!!!!!

        int cnt = 6 * (r - l).abs() / (translation1.abs() + translation2.abs() + translation3.abs());
        for (int nx = -cnt; nx <= cnt; nx++)
            for (int ny = -cnt; ny <= cnt; ny++)
                for (int nz = -cnt; nz <= cnt; nz++)
                {

                    Vector pb = translation1 * nx + translation2 * ny + translation3 * nz;
                    for (auto &vector_type : latt)
                    {

                        auto &[vec, a] = vector_type;

                        Vector pa = pb + translations.to_euclidus(vec);

                        if (!geometry::IsVectorInCube(pa, p0, p1))
                            continue;

                        insert(pa, std::make_shared<grid::atom::Atom>(a));
                    }
                }
    };

    void Grid::ClearParallelep(const geometry::Vector &l, const geometry::Vector &r)
    {
        for (double x = l.x; x < r.x; x += size_chunk)
            for (double y = l.y; y < r.y; y += size_chunk)
                for (double z = l.z; z < r.z; z += size_chunk)
                {
                    auto p = calcChunkPos(Vector(x, y, z));
                    auto pvec = chunks.find(p);
                    if (pvec == chunks.end())
                        continue;
                    auto &chunk = *pvec->second;
                    for (auto a = chunk.begin();a!=chunk.end();)
                    {
                        auto& pos = a->first;
                        a++;
                        if (geometry::IsVectorInCube(pos, l, r))
                        {

                            chunk.erase(pos);
                        }
                    }
                }
    }

    double Grid::Chunk_Size() const { return size_chunk; };

    size_t Grid::count() const
    
    {
    
        size_t cnt = 0;
    
        for (auto a : chunks)
    
        {
    
            cnt += a.second->count();
    
        };
    
        return cnt;
    
    };

    Vector Grid::getCycledVector(const Vector &v) const
    {
        if (geometry::IsVectorInCube(v, llim, rlim))
            return v;
        auto [x, y, z] = v - llim;
        auto [px, py, pz] = rlim - llim;
        if (cycle_x)
            x = math::modulo(x, px);
        if (cycle_y)
            y = math::modulo(y, py);
        if (cycle_z)
            z = math::modulo(z, pz);
        return {x, y, z};
    }

    string Grid::to_xyz(double mult) const
    {
        std::stringstream ss(std::ios_base::out);
        ss << count() << "\n\n";

        for_each([&ss,mult](const Vector& pos, const std::shared_ptr<atom::Atom>& atom)mutable
        {
            Vector multiplied = pos * mult;
            const auto &material = atom->Material();
            ss << material->Name() << " " << multiplied.x << " " << multiplied.y << " " << multiplied.z << "\n";
        });

        return ss.str();
    }

    string Grid::to_xyz() const
    {
        return to_xyz(1);
    }

    void Grid::from_xyz(const string &s, double div)
    {
        std::stringstream ss(s, std::ios_base::in);
        size_t count;
        ss >> count;
        auto tend = types.end();
        for (size_t i = 0; i < count; i++)
        {
            string name;
            double x, y, z;
            ss >> name >> x >> y >> z;
            auto type = types.find(name);
            assert(type != tend);
            auto &T = type->second;
            Vector pos(x, y, z);
            pos /= div;
            insert(pos, std::make_shared<atom::Atom>(T));
            if (rlim.x < x)
            {
                rlim.x = x;
            }
            if (rlim.y < y)
            {
                rlim.y = y;
            }
            if (rlim.z < z)
            {
                rlim.z = z;
            }
        }
    }

    void Grid::from_xyz(const string &s)
    {
        from_xyz(s, 1);
    }

    
    void Grid::for_each(for_each_const_callbak callback) const
    {
        for(auto& [pos,chunk]:chunks) 
        {
            chunk->for_each(callback);
        }
    };

    void Grid::for_each(const Vector& pos, double dist, for_each_const_callbak callback) const
    {
        int translation_x, translation_y, translation_z;
        int mtranslation_x, mtranslation_y, mtranslation_z;

        Vector translations = Sizes();

        size_t chunksAround;

        double x, y, z;
        double mx, my, mz;

        Vector c0 = calcChunkPos(pos);

        Vector cllim = c0;
        Vector crlim = c0 + Vector(size_chunk, size_chunk, size_chunk);

        Vector rp = pos + Vector(dist, dist, dist);
        Vector lp = pos - Vector(dist, dist, dist);

        mtranslation_x = 0;
        mtranslation_y = 0;
        mtranslation_z = 0;

        chunksAround = std::floor(dist / size_chunk) + 1;

        mx = chunksAround * size_chunk;
        my = chunksAround * size_chunk;
        mz = chunksAround * size_chunk;

        if (geometry::IsVectorInCube(rp, cllim, crlim) && geometry::IsVectorInCube(lp, cllim, crlim))
        {
            mx = 0;
            my = 0;
            mz = 0;
        }

        if (Cyclic<'x'>())
        {
            mtranslation_x = 1;
        }
        if (Cyclic<'y'>())
        {
            mtranslation_y = 1;
        }
        if (Cyclic<'z'>())
        {
            mtranslation_z = 1;
        }


        for (x=-mx; x <= mx; x += size_chunk)
        {
            for (y=-my; y <= my; y += size_chunk)
            {
                for (z=-mz; z <= mz; z += size_chunk)
                {
                    //Вектор переноса к следующему чанку.
                    Vector delta(x,y,z); 


                    for (translation_x = -mtranslation_x; translation_x <= mtranslation_x; translation_x++)
                    {
                        for (translation_y = -mtranslation_y; translation_y <= mtranslation_y; translation_y++)
                        {
                            for (translation_z = -mtranslation_z; translation_z <= mtranslation_z; translation_z++)
                            {
                                //Вектор, для выбора чанка с учетом периодических границ и трансляций.
                                Vector PeriodicTranslation = translations.coord_mul({(double)(translation_x), (double)(translation_y), (double)(translation_z)});
                                Vector ChunkTranslation = delta+PeriodicTranslation;
                                //Проверка, что хоть одна точка чанка достаточно близка к центру сферы поиска
                                if ((ChunkTranslation).abs()-size_chunk < dist && getChunkSimple(pos+delta)) 
                                {
                                    Vector ChunkTranslated = pos + ChunkTranslation;

                                    auto pChunk = getChunk(ChunkTranslated);
                                    if(pChunk) {
                                        pChunk->for_each([&](const auto& PosAtom, auto& atom) mutable
                                        {
                                            auto deltaPos = PosAtom-pos+PeriodicTranslation;
                                            double distanceTranslated = deltaPos.abs();
                                            if(distanceTranslated < dist) 
                                            {
                                                callback(PosAtom,atom);
                                            }
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}