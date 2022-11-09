//
// Created by Danil-Ponkratov
//

#define _USE_MATH_DEFINES
#include <params.hpp>

#include <settings.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>

#include <sgs.hpp>

params::params()
    : Ea1(settings["Ea1"]),
    Ea2(settings["Ea2"]),
    Ea3(settings["Ea3"]),
    Ea4(settings["Ea4"]),
    ER1(settings["ER1"]),
    ER2(settings["ER2"]),
    ER3(settings["ER3"]),
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
    Temperature(settings["Temperature"]),
    vac_size(settings["vac_size"])
{
  std::cout<<"params initializing"<<std::endl;
  sina = sin(M_PI*alpha/180);
  sinb = sin(M_PI*beta/180);
  sinc = sin(M_PI*gamma/180);
  cosa = sin(M_PI*alpha/180);
  cosb = sin(M_PI*beta/180);
  cosc = sin(M_PI*gamma/180);
  Ea1=Ea1*ELVOLT;
  Ea2=Ea2*ELVOLT;
  Ea3=Ea3*ELVOLT;
  Ea4=Ea4*ELVOLT;
  ER1=ER1*ELVOLT;
  ER2=ER2*ELVOLT;
  ER3=ER3*ELVOLT;
}

params::params(double _E)
    : Ea1(settings["Ea1"]),
    Ea2(settings["Ea2"]),
    Ea3(settings["Ea3"]),
    Ea4(settings["Ea4"]),
    ER1(settings["ER1"]),
    ER2(settings["ER2"]),
    ER3(settings["ER3"]),
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
    Temperature(settings["Temperature"]),
    vac_size(settings["vac_size"])
{
  std::cout<<"params initializing"<<std::endl;
  sina = sin(M_PI*alpha/180);
  sinb = sin(M_PI*beta/180);
  sinc = sin(M_PI*gamma/180);
  cosa = sin(M_PI*alpha/180);
  cosb = sin(M_PI*beta/180);
  cosc = sin(M_PI*gamma/180);

  Ea1=Ea1*ELVOLT;
  Ea2=Ea2*ELVOLT;
  Ea3=Ea3*ELVOLT;
  Ea4=Ea4*ELVOLT;
  ER1=ER1*ELVOLT;
  ER2=ER2*ELVOLT;
  ER3=ER3*ELVOLT;
}