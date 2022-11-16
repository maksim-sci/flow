#include <geometry/geometry.hpp>
namespace geometry
{
    Vector Geometry::to_euclidus(const Vector &rhs) const
    {

        return (a * rhs.x + b * rhs.y + c * rhs.z);
    }

    Vector Geometry::from_euclidus(const Vector &rhs) const
    {

        Vector a_t(a.x, b.x, c.x);

        Vector b_t(a.y, b.y, c.y);

        Vector c_t(a.z, b.z, c.z);

        double av = a_t.norm().scalar_mul(rhs) / a_t.abs();

        double bv = b_t.norm().scalar_mul(rhs) / b_t.abs();

        double cv = c_t.norm().scalar_mul(rhs) / c_t.abs();

        return Vector(av, bv, cv);
    }
}