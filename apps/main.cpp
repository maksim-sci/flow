//
// Created by Danil-Ponkratov
//

#include <app/app.hpp>

#include <tuple>

int main() {
  auto limGrid = std::make_tuple(15, 10, 20);
  auto limLug = std::make_tuple(5, 7, 3, 6, 0, 20);
  bool generateGrid = true; bool isFromFile = false; bool isToFile = true;
  int cntFlow = 1000; int freqRecording = 50; int msgDebug = 50;
  auto r = runFlow(limGrid, limLug, generateGrid, isFromFile, isToFile);
  r.run(cntFlow, freqRecording, msgDebug);
  return 0;
}
