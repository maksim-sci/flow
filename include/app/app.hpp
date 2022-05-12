//
// Created by Danil-Ponkratov
//

#pragma once

#include <grid/grid.hpp>
#include <section/section.hpp>
#include <params/params.hpp>

#include <string>
#include <tuple>

class runFlow
{
public:
  runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
          bool generateGrid, bool isFromFile, bool isToFile);

  void run(int cntFlow, int freqRecording, int msgDebug);

private:
  class grid grid;
  section sec;
  params param;
  bool isToFile;
};
