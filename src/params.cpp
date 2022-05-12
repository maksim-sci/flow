//
// Created by Danil-Ponkratov
//

#include <params/params.hpp>

params::params()
    : Ea1(1.0),
      Ea2(0.9),
      Ea3(0.82),
      Ea4(0.82),
      y(7.0),
      a(3.589),
      T(300.0),
      k(8.617e-5),
      E(5.9e6),
      e(4.803e-10),
      v(10.0e12)
{}

params::params(double E)
    : Ea1(1.0),
      Ea2(0.9),
      Ea3(0.82),
      Ea4(0.82),
      y(7.0),
      a(3.589),
      T(300.0),
      k(8.617e-5),
      e(4.803e-10),
      v(10.0e12),
      E(E)
{}