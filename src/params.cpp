//
// Created by Danil-Ponkratov
//

#include <params/params.hpp>

#include <../include/json.hpp>
#include <string>
#include <fstream>
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
{}

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