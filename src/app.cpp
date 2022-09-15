//
// Created by Danil-Ponkratov
//
//
//Upgraded by Maksim Solovev
//
#include "types.h"
#include <app/app.hpp>
#include <flow/flow.hpp>
#include <memory>
#include <sstream>
#include <locale>
#include <iomanip>
#include <filesystem>
#include <istream>
#include <fstream>

namespace fs = std::filesystem;
runFlow::runFlow(std::tuple<int, int, int> limGrid, std::tuple<int, int, int, int, int, int> limLug,
          bool generateGrid, std::string from_file, std::string _outperiodic,bool _uniqe_folder)
    : grid(std::get<0>(limGrid), std::get<1>(limGrid), std::get<2>(limGrid), &param,generateGrid)
      , sec(section())
      ,outperiodic(_outperiodic),
      param(),
      f(grid, sec, param),
      step_number(0),
      freqRecording(999999999)
{
  if (from_file!="")
    fromFile(from_file, grid);
  if(outperiodic!="" && _uniqe_folder)
  {

    std::time_t time = std::time(nullptr);
    {
      using std::to_string;
      using std::get;
      char mbstr[100];
      std::strftime(mbstr, sizeof(mbstr), "%H_%M_%S(%d_%m_%Y)", std::localtime(&time));
      outperiodic = (fs::absolute(fs::path(outperiodic))
      .append(to_string(get<0>(limGrid))+"_"+
      to_string(get<1>(limGrid))+"_"+
      to_string(get<2>(limGrid))+"_"+
      std::string(mbstr)
      )).string();
    }
    fs::create_directories(outperiodic);
    auto gridout_path = fs::path(outperiodic).append("grid");
    gridOut = gridout_path.string();

    fs::create_directories(gridout_path);
    fs::copy_file("settings.json",fs::path(outperiodic).append("settings.json"));

  }
  auto [lx, rx, ly, ry, lz, rz] = limLug;
  grid.addLug(lx, rx, ly, ry, lz, rz);
}
void runFlow::step()
{
  if ( (step_number % freqRecording == 0) )
  {
    if(outperiodic!="")
    {
      std::ofstream outf(fs::absolute(outperiodic).append("counts"),std::ios_base::app);
      if(outf)
      {
        for(int i = 0;i<5;i++)
        {
          outf<<(TypeReaction)i<<":"<<f.statistic[i]<<" ";
          f.statistic[i]=0;
        }
        outf<<std::endl;
        outf.close();
      }
      std::cout << "recording in flow " << step_number << std::endl;
      toFile(fs::absolute(gridOut).append("run_"+std::to_string(step_number)+".xyz").string(), grid);
      std::ofstream current_stream(fs::absolute(outperiodic).append("current.txt"),std::ios_base::app);
      current_stream<<step_number<<"\t"<<f.I<<std::endl;
      f.I = 0;
      current_stream.close();
    }
  }
  f.run();

  step_number++;
  sec.clear();
}

void runFlow::run(int cntFlow)
{
  for (int i = 0; i < cntFlow; i++)
  {
    step();
  }

  if (cntFlow < freqRecording)

  f.getStatistic();
}