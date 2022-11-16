#include <grid/lattice.hpp>

namespace grid {
    void Lattice::add(Vector pos, std::shared_ptr<Type> atype)
    
    {
    
        if (positions.find(pos) != positions.end())
    
        {
    
            throw fmt::format("{}:{} position: {} is already occuped", __FILE__, __LINE__, pos);
    
        }
    
        types.insert(atype);
    
        positions.insert({pos, atype});
    
    };
}