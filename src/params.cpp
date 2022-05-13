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
  // std::cout<<settings["Ea1"]<<std::endl;
  //   std::cout<<settings["Ea2"]<<std::endl;
  //   std::cout<<settings["Ea3"]<<std::endl;
  //   std::cout<<settings["Ea4"]<<std::endl;
  //   std::cout<<settings["y"]<<std::endl;
  //   std::cout<<settings["a"]<<std::endl;
  //   std::cout<<settings["T"]<<std::endl;
  //   std::cout<<settings["k"]<<std::endl;
  //   std::cout<<settings["E"]<<std::endl;
  //   std::cout<<settings["e"]<<std::endl;
  //   std::cout<<settings["v"]<<std::endl;

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