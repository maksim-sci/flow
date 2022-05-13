//
// Created by Danil-Ponkratov
//

#include <params/params.hpp>

#include <../include/json.hpp>
#include <string>
#include <fstream>
#include <iostream>
auto settings = nlohmann::json::parse(std::ifstream("settings.json"));

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
    v(settings["v"])
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
    E(_E)
{}