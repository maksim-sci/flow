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
        for (auto a = begin(); !a.Finished();)
        {
            auto& atom = a.aiter->second;
            auto pos = a.aiter->first;
            pos *= mult;
            const auto &material = atom->Material();
            ss << material->Name() << " " << pos.x << " " << pos.y << " " << pos.z << "\n";
            ++a;
        };

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

    Grid::GridIteratorOnceSinglePass::GridIteratorOnceSinglePass(const Grid &g)
    {
        citer = g.chunks.begin();
        citerend = g.chunks.end();
        aiter = citer->second->begin();
        aiterend = citer->second->end();
    };
    
    Grid::GridIteratorOnceSinglePass &Grid::GridIteratorOnceSinglePass::operator++()
    {
        if (aiter != aiterend)
        {
            aiter++;
        }
        while(aiter==aiterend && (++citer)!=citerend) {
            aiter = citer->second->begin();
        }
        return *this;
    };

    void Grid::GridIteratorOnceSinglePass::end()
    
    {
    
        citer = citerend;
    
        aiter = aiterend;
    
    };

    void Grid::GridIteratorDistLimSinglePass::iterateTillCorrectOrEnd()
    {
        // what the hell?
        //Да, это самая уродливая функция, которую вы видели, мне тоже очень нравится.
        for (; x <= mx; x += chunk_size)
        {
            for (; y <= my; y += chunk_size)
            {
                for (; z <= mz; z += chunk_size)
                {
                    Vector delta(x,y,z);

                    for (; translation_x <= mtranslation_x; translation_x++)
                    {
                        for (; translation_y <= mtranslation_y; translation_y++)
                        {
                            for (; translation_z <= mtranslation_z; translation_z++)
                            {
                                //fmt::print("x,y,z: ({} {} {})  tx,ty,tz: ({} {} {})\n",x,y,z,translation_x,translation_y,translation_z);
                                if (delta.abs()-chunk_size < distance)
                                {
                                    Vector ptrans = pos + translations.coord_mul(Vector((double)(translation_x), (double)(translation_y), (double)(translation_z)));
                                    Vector p1 = ptrans+ delta;

                                    Vector p0 = grid->calcChunkPos(p1);
                                    auto chunk = chunks->find(p0);
                                    if (chunk != chunks->end())
                                    {

                                        aiter = chunk->second->begin();
                                        aiterend = chunk->second->end();
                                        while(aiter!=aiterend) {
                                            auto& Vectorp = aiter->first;
                                            if((Vectorp - ptrans).abs() < distance) {
                                                return;
                                            }
                                            aiter++;
                                        }
                                    }
                                }
                            }
                            translation_z = -mtranslation_z;
                        }
                        translation_y = -mtranslation_y;
                    }
                    translation_x = -mtranslation_x;
                }
                z = pos.z - chunksAround * chunk_size;
            }
            y = pos.y - chunksAround * chunk_size;
        }
        finished = true;
    }

    Grid::GridIteratorDistLimSinglePass::GridIteratorDistLimSinglePass(Grid &g, double dist, const Vector &p) : chunks(&g.chunks), finished(false), pos(p)
    {
        // check 1 chunk
        c0 = g.calcChunkPos(p);

        chunk_size = g.Chunk_Size();

        chunks = &g.chunks;

        grid = &g;

        distance = dist;

        Vector cllim = c0;
        Vector crlim = c0 + Vector(chunk_size, chunk_size, chunk_size);

        Vector rp = p + Vector(dist, dist, dist);
        Vector lp = p - Vector(dist, dist, dist);

        translation_x = 0;
        mtranslation_x = 0;
        translation_y = 0;
        mtranslation_y = 0;
        translation_z = 0;
        mtranslation_z = 0;

        if (geometry::IsVectorInCube(rp, cllim, crlim) && geometry::IsVectorInCube(lp, cllim, crlim))
        {
            mx = 0;
            x = 0;
            y = 0;
            my = 0;
            z = 0;
            mz = 0;
            auto chunk = grid->chunks.find(c0);

            aiter = chunk->second->begin();
            aiterend = chunk->second->end();
            if (aiter == aiterend)
                finished = true;
            return;
        }
        chunksAround = std::floor(dist / chunk_size) + 1;
        x = pos.x - chunksAround * chunk_size;
        z = pos.z - chunksAround * chunk_size;
        y = pos.y - chunksAround * chunk_size;
        mx = pos.x + chunksAround * chunk_size;
        my = pos.y + chunksAround * chunk_size;
        mz = pos.z + chunksAround * chunk_size;

        if (g.Cyclic<'x'>())
        {
            translation_x = -1;
            mtranslation_x = 1;
        }
        if (g.Cyclic<'y'>())
        {
            translation_y = -1;
            mtranslation_y = 1;
        }
        if (g.Cyclic<'z'>())
        {
            translation_z = -1;
            mtranslation_z = 1;
        }

        translations = g.Sizes();

        iterateTillCorrectOrEnd();
    };
    bool Grid::GridIteratorDistLimSinglePass::operator==(const Grid::GridIteratorDistLimSinglePass &g) const
    {
        if (x == mx && y == my && z == mz)
            if (aiter == aiterend)
                return true;
        return false;
    };

    Grid::GridIteratorDistLimSinglePass &Grid::GridIteratorDistLimSinglePass::operator++()
    {
        aiter++;
        if (aiter == aiterend)
        {
            //cringe
            if(translation_x==mtranslation_x) {
                if(translation_y==mtranslation_y) 
                    if(translation_z==mtranslation_z)
                        if(z==mz)
                            if(y==my)
                                if(x==mx) {
                                    finished = true;
                                    return *this;
                                }
                                else x+=chunk_size;
                            else y+=chunk_size;
                        else z+=chunk_size;
                    else translation_z++;
                else translation_y++;
            }
            else {
                translation_x++;
            }
            iterateTillCorrectOrEnd();
        }
        return *this;
    };

}