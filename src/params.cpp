//
// Created by Danil-Ponkratov
//

#include <params/params.hpp>

#include <settings.hpp>
#include <string>
#include <fstream>
#include <iostream>

params::params()
    : Ea1(settings["Ea1"]),
    Ea2(settings["Ea2"]),
    Ea3(settings["Ea3"]),
    Ea4(settings["Ea4"]),
    y(settings["y"]),
    a(settings["a"]),
    T(settings["T"]),
    k(settings["k"]),
    E(settings["E"]),
    e(settings["e"]),
    v(settings["v"]),
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
  // std::cout<<Ea1<<std::endl;
  // std::cout<<Ea2<<std::endl;
  // std::cout<<Ea3<<std::endl;
  // std::cout<<Ea4<<std::endl;
  // std::cout<<y<<std::endl;
  // std::cout<<a<<std::endl;
  // std::cout<<T<<std::endl;
  // std::cout<<k<<std::endl;
  // std::cout<<E<<std::endl;
  // std::cout<<e<<std::endl;
  // std::cout<<v<<std::endl;

}

params::params(double _E)
    : Ea1(settings["Ea1"]),
    Ea2(settings["Ea2"]),
    Ea3(settings["Ea3"]),
    Ea4(settings["Ea4"]),
    y(settings["y"]),
    a(settings["a"]),
    T(settings["T"]),
    k(settings["k"]),
    e(settings["e"]),
    v(settings["v"]),
    E(_E),
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
{}