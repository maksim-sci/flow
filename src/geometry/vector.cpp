#include <math.h>
#include <fmt/format.h>
#include "vector.hpp"

namespace geometry
{
  bool IsVectorInCube(const Vector &a, const Vector &l, const Vector &r)
  {
    if ((l.x > a.x) || (l.y > a.y) || (l.z > a.z))
      return false;
    if ((r.x < a.x) || (r.y < a.y) || (r.z < a.z))
      return false;

    return true;
  }

  Vector &Vector::operator=(const Vector &rhs)
  {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    return *this;
  };
  Vector &Vector::operator-=(const Vector &rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  };
  Vector &Vector::operator+=(const Vector &rhs)
  {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  };
  Vector &Vector::operator*=(const double mul)
  {
    x *= mul;
    y *= mul;
    z *= mul;
    return *this;
  };
  Vector &Vector::operator/=(const double mul)
  {
    x /= mul;
    y /= mul;
    z /= mul;
    return *this;
  };
  Vector Vector::operator+(const Vector &rhs) const
  {
    Vector A = *this;
    return A += rhs;
  };
  Vector Vector::operator-(const Vector &rhs) const
  {
    Vector A = *this;
    return A -= rhs;
  };
  Vector Vector::operator*(const double mul) const
  {
    Vector A = *this;
    return A *= mul;
  };
  Vector Vector::operator/(const double mul) const
  {
    Vector A = *this;
    return A /= mul;
  };
  bool Vector::operator<(const Vector &rhs) const
  {
    if (x < rhs.x)
      return true;
    if (x > rhs.x)
      return false;
    if (y < rhs.y)
      return true;
    if (y > rhs.y)
      return false;
    if (z < rhs.z)
      return true;
    return false;
  }
  bool Vector::operator>(const Vector &rhs) const
  {
    if (x > rhs.x)
      return true;
    if (x < rhs.x)
      return false;
    if (y > rhs.y)
      return true;
    if (y < rhs.y)
      return false;
    if (z > rhs.z)
      return true;
    return false;
  }
  Vector Vector::norm() const { return *this / abs(); };
  double Vector::abs() const { return powf(x * x + y * y + z * z, 0.5); };
  double Vector::scalar_mul(const Vector &rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; };
  Vector Vector::coord_mul(const Vector &rhs) const { return {x * rhs.x, y * rhs.y, z * rhs.z}; };
  bool Vector::operator==(const Vector &rhs) const
  {
    auto delta = *this - rhs;
    double modd = delta.abs();
    double mod1 = this->abs();
    double mod2 = rhs.abs();
    return modd <= (fmin(mod1, mod2)) * 1e-6;
  };
  bool Vector::operator!=(const Vector &rhs) const { return !(rhs == *this); };

}