//
// Created by Danil-Ponkratov
//

#include <app/app.hpp>
#include <flow/flow.hpp>

std::string out("C:/Users/user/Desktop/New Albany/run_");
runFlow::runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
                 bool generateGrid, bool isFromFile, bool isToFile)
    : grid(std::get<0>(limGrid), std::get<1>(limGrid), std::get<2>(limGrid), generateGrid)
      , sec(section())
      , param(params())
      , isToFile(isToFile)
{
  if (isFromFile)
    fromFile("data/after_run.xyz", grid);

  auto [lx, rx, ly, ry, lz, rz] = limLug;
  grid.addLug(lx, rx, ly, ry, lz, rz);
}


void runFlow::run(int cntFlow, int freqRecording, int msgDebug)
{
  auto f = flow(grid, sec, param);
  for (int i = 0; i < cntFlow; i++)
  {
    f.run();
    if ( isToFile && (i % freqRecording == 0) )
    {
      std::cout << "recording in flow " << i << std::endl;
      toFile((out+std::to_string(i)+".xyz"), grid);
      toFile("data/after_run.xyz", grid);
    }
    sec.clear();

    if ( msgDebug != 0 && (i % msgDebug == 0) )
      std::cout << "flow " << i << std::endl;
  }

  if (isToFile && cntFlow < freqRecording)
    toFile("data/after_run.xyz", grid);

  f.getStatistic();
}