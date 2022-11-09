//
// Created by Danil-Ponkratov
//

#pragma once

#include "types.h"
#include <grid.hpp>
#include <section.hpp>
#include <params.hpp>
#include <flow.hpp>

#include <string>
#include <tuple>
#include <memory>

class runFlow
{
public:
  runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
          bool generateGrid, std::string from_file, std::string _outperiodic,bool _unique_folder);

  void step();
  void run(int cntFlow);
  params param;

  size_t freqRecording;
  size_t step_number;
  std::string outperiodic;
  std::string gridOut;
  class grid grid;
  section sec;
  bool isToFile;
  flow f;
private:
};
