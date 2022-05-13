//
// Created by Danil-Ponkratov
//

#pragma once

#include <grid/grid.hpp>
#include <section/section.hpp>
#include <params/params.hpp>

#include <string>
#include <tuple>
#include <memory>

class runFlow
{
public:
  runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
          bool generateGrid, std::string from_file, std::string _outperiodic,bool _unique_folder);

  void run(int cntFlow, int freqRecording, int msgDebug);

  

private:
  std::string outperiodic;
  class grid grid;
  section sec;
  params param;
  bool isToFile;
};
