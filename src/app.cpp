//
// Created by Danil-Ponkratov
//
//
//Upgraded by Maksim Solovev
//

#include <app/app.hpp>
#include <flow/flow.hpp>
#include <memory>
runFlow::runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
          bool generateGrid, std::string from_file, std::string _outperiodic)
    : grid(std::get<0>(limGrid), std::get<1>(limGrid), std::get<2>(limGrid), generateGrid)
      , sec(section())
      ,outperiodic(_outperiodic),
      param()
{
  if (from_file!="")
    fromFile(from_file, grid);

  auto [lx, rx, ly, ry, lz, rz] = limLug;
  grid.addLug(lx, rx, ly, ry, lz, rz);
}


void runFlow::run(int cntFlow, int freqRecording, int msgDebug)
{
  auto f = flow(grid, sec, param);
  for (int i = 0; i < cntFlow; i++)
  {
    f.run();
    if ( (i % freqRecording == 0) )
    {
      if(outperiodic!="")
        std::cout << "recording in flow " << i << std::endl;
        toFile((outperiodic+std::to_string(i)+".xyz"), grid);
    }
    sec.clear();
  }

  if (cntFlow < freqRecording)

  f.getStatistic();
}