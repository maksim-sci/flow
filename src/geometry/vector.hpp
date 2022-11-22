#pragma once
#include <math.h>
#include <fmt/format.h>
namespace geometry {
    class Vector {
    public:
        double x;
        double y;
        double z;
    public:
        inline Vector():x(0),y(0),z(0){};
        inline Vector(double _x, double _y, double _z):x(_x),y(_y),z(_z) {};

        inline ~Vector() {};

        Vector& operator=(const Vector& rhs);
        Vector& operator-=(const Vector& rhs);
        Vector& operator+=(const Vector& rhs);
        Vector& operator*=(const double mul);
        Vector& operator/=(const double mul);
        Vector operator+ (const Vector& rhs) const;
        Vector operator- (const Vector& rhs) const;
        Vector operator* (const double mul) const;
        Vector operator/ (const double mul) const;
        bool operator<(const Vector& rhs) const;
        bool operator>(const Vector& rhs) const;
        Vector norm() const;
        double abs()  const;
        double scalar_mul (const Vector& rhs) const;
        Vector coord_mul (const Vector& rhs) const;
        bool operator== (const Vector& rhs) const;
        bool operator!= (const Vector& rhs) const;


    };

  bool IsVectorInCube(const Vector& a,const  Vector& l,const  Vector& r);
    
}

template <> struct fmt::formatter<geometry::Vector> {
  char presentation = 'f';
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'f' || *it == 'e')) presentation = *it++;
    if (it != end && *it != '}') throw format_error("invalid format");
    return it;
  }
  template <typename FormatContext>
  inline auto format(const geometry::Vector& p, FormatContext& ctx) const -> decltype(ctx.out()) {
    return presentation == 'f'
              ? fmt::format_to(ctx.out(), "Vector({:f}, {:f}, {:f})", p.x, p.y,p.z)
              : fmt::format_to(ctx.out(), "Vector({:.4e}, {:.4e}, {:.4e})", p.x, p.y,p.z);
  }
};

namespace std {
  template <> struct hash<geometry::Vector> {
    size_t operator()(const geometry::Vector& vec) const {
        auto a = [](double d) {return std::hash<double>()(d);};
        return a(vec.x)^a(vec.y)^a(vec.z);
    }
  };
}