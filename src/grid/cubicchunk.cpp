#include "cubicchunk.hpp"
#include <algorithm>

using geometry::Vector;


namespace grid {
    namespace chunk {
        CubicChunk::CubicChunk() : pos(0, 0, 0), size(0), atoms(){};

        CubicChunk::CubicChunk(double x, double y, double z, double size) : pos(x, y, z), size(size), atoms(){};

        CubicChunk::CubicChunk(const Vector &pos, double size) : pos(pos), size(size), atoms(){};

        size_t CubicChunk::cnt() const { return atoms.size(); };

        Vector CubicChunk::getLeftLimit() const { return pos; };

        Vector CubicChunk::getRightLimit() const { return (pos + Vector(size, size, size)); };

        double CubicChunk::dim_size() const { return size; };

        bool CubicChunk::inPosIn(const Vector &v) const
        {
            const auto &l = pos;
            const auto r = getRightLimit();
            return geometry::IsVectorInCube(v, l, r);
        };

        size_t CubicChunk::count() const
        {
            size_t cnt = 0;
            for (auto a : atoms)
            {
                cnt++;
            };
            return cnt;
        };

        std::shared_ptr<atom::Atom> CubicChunk::get(const Vector &pos) const
        
        {
        
            auto result = atoms.find(pos);
        
            if (result != atoms.end())
        
            {
        
                return result->second;
        
            }
        
            return nullptr;
        
        };

        InsertionResults CubicChunk::insert(Vector pos, std::shared_ptr<grid::atom::Atom> atom) {
            auto a = atoms.find(pos);
            if(a!=atoms.end()){return InsertionResults::RepeatedInsertion;}
            if(!inPosIn(pos)) {return InsertionResults::InsertionOutOfBounds;}
            atoms.insert({pos,atom}); 
            return InsertionResults::OK; 
        };

        bool CubicChunk::erase(const Vector &pos)
        
        {
        
            atoms.erase(pos);
        
            return true;
        
        };

        void CubicChunk::for_each(std::function<void(const Vector&,std::shared_ptr<atom::Atom>&)> callback)
        {
            for(auto data:atoms) {
                auto& [pos,atom] = data;
                callback(pos,atom);
            }
        };
    }
}
