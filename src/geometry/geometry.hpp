#pragma once
#include "vector.hpp"

namespace geometry {
    class Geometry {
        Vector a;
        Vector b;
        Vector c;

    public: 


        inline Geometry():a(1.,1.,1.),b(1.,1.,1.),c(1.,1.,1.) { };
        inline Geometry(const Vector& _a,const Vector& _b,const Vector& _c):a(_a),b(_b),c(_c) {};

        inline ~Geometry(){};

        Vector to_euclidus(const Vector& rhs) const;

        Vector from_euclidus(const Vector&rhs) const;

        inline const Vector& A() const {return a;};
        inline const Vector& B() const {return b;};
        inline const Vector& C() const {return c;};
    };

}