#pragma once
#include <../include/json.hpp>
#include <fstream>
#include <settings.hpp>
nlohmann::json settings = nlohmann::json::parse(std::ifstream("settings.json"));
