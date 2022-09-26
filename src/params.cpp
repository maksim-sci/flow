//
// Created by Danil-Ponkratov
//

#define _USE_MATH_DEFINES
#include <params/params.hpp>

#include <settings.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>

params::params()
    : Ea1(settings["Ea1"]),
    Ea2(settings["Ea2"]),
    Ea3(settings["Ea3"]),
    Ea4(settings["Ea4"]),
    y(settings["y"]),
    a(settings["a"]),
    b(settings["b"]),
    c(settings["c"]),
    alpha(settings["alpha"]),
    beta(settings["beta"]),
    gamma(settings["gamma"]),
    T(settings["T"]),
    k(settings["k"]),
    E(settings["E"]),
    e(settings["e"]),
    v(settings["v"]),
    U(settings["U"]),
    R1mult(settings["R1mult"]),
    R2mult(settings["R2mult"]),
    R3mult(settings["R3mult"]),
    R4mult(settings["R4mult"]),
    E1mult(settings["E1mult"]),
    E2mult(settings["E2mult"]),
    E3mult(settings["E3mult"]),
    AE1(settings["AE1"]),
    AE2(settings["AE2"]),
    AE3(settings["AE3"]),
    DE2(settings["DE2"]),
    DE3(settings["DE3"]),
    le(settings["le"]),
    hconst(settings["hconst"]),
    kb(settings["kb"]),
    Temperature(settings["Temperature"])
{

  sina = sin(M_PI*alpha/180);
  sinb = sin(M_PI*beta/180);
  sinc = sin(M_PI*gamma/180);
  cosa = sin(M_PI*alpha/180);
  cosb = sin(M_PI*beta/180);
  cosc = sin(M_PI*gamma/180);

}

params::params(double _E)
    : Ea1(settings["Ea1"]),
    Ea2(settings["Ea2"]),
    Ea3(settings["Ea3"]),
    Ea4(settings["Ea4"]),
    y(settings["y"]),
    a(settings["a"]),
    b(settings["b"]),
    c(settings["c"]),
    alpha(settings["alpha"]),
    beta(settings["beta"]),
    gamma(settings["gamma"]),
    T(settings["T"]),
    k(settings["k"]),
    e(settings["e"]),
    v(settings["v"]),
    E(_E),
    U(settings["U"]),
    R1mult(settings["R1mult"]),
    R2mult(settings["R2mult"]),
    R3mult(settings["R3mult"]),
    R4mult(settings["R4mult"]),
    E1mult(settings["E1mult"]),
    E2mult(settings["E2mult"]),
    E3mult(settings["E3mult"]),
    AE1(settings["AE1"]),
    AE2(settings["AE2"]),
    AE3(settings["AE3"]),
    DE2(settings["DE2"]),
    DE3(settings["DE3"]),
    le(settings["le"]),
    hconst(settings["hconst"]),
    kb(settings["kb"]),
    Temperature(settings["Temperature"])
{
  sina = sin(M_PI*alpha/180);
  sinb = sin(M_PI*beta/180);
  sinc = sin(M_PI*gamma/180);
  cosa = sin(M_PI*alpha/180);
  cosb = sin(M_PI*beta/180);
  cosc = sin(M_PI*gamma/180);
}