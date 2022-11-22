#pragma once
#include <cmath>

namespace math {
    template <typename T>
    T modulo(const T& a,const  T& b) {
        return a-std::floor(a/b)*b;
    }
}