//
// Created by Danil-Ponkratov
//
//
//Upgraded by Maksim Solovev
//

#include <app/app.hpp>
#include <flow/flow.hpp>
#include <memory>
#include <sstream>
#include <locale>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;
runFlow::runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
          bool generateGrid, std::string from_file, std::string _outperiodic,bool _uniqe_folder)
    : grid(std::get<0>(limGrid), std::get<1>(limGrid), std::get<2>(limGrid), generateGrid)
      , sec(section())
      ,outperiodic(_outperiodic),
      param()
{
  if (from_file!="")
    fromFile(from_file, grid);
  if(outperiodic!="" && _uniqe_folder)
  {

    std::time_t time = std::time(nullptr);
    outperiodic = fs::absolute(fs::absolute(fs::path(outperiodic)).append(std::to_string(time))).string();
    fs::create_directories(outperiodic);
    fs::copy_file("settings.json",fs::path(outperiodic).append("settings.json"));

  }
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
        toFile(fs::absolute(outperiodic).append("run_"+std::to_string(i)+".xyz").string(), grid);
    }
    sec.clear();

    //f.reactionsBox.clear();
  }

  if (cntFlow < freqRecording)

  f.getStatistic();
}