//
// Created by Danil-Ponkratov
//

#include <app/app.hpp>


#include <fstream>
#include <string>
#include <iostream>

#include <../include/json.hpp>

std::string path_to_settings = "settings.json";
using namespace std;

namespace settings
{
    int sizex,sizey,sizez;
    int minx,miny,minz,maxx,maxy,maxz;
    int recordPeriod,timeLimit;
    bool saveToFile,generateGrid,unique_folder;
    std::string from_file,out_folder;

}


void loadsettings()
{
    std::ifstream ifs(path_to_settings);
    auto settings = nlohmann::json::parse(ifs);
    ifs.close();
    settings::sizex = settings["sizex"];
    settings::sizey = settings["sizey"];
    settings::sizez = settings["sizez"];
    settings::minx = settings["minx"];
    settings::miny = settings["miny"];
    settings::minz = settings["minz"];
    settings::maxx = settings["maxx"];
    settings::maxy = settings["maxy"];
    settings::maxz = settings["maxz"];
    settings::recordPeriod = settings["recordPeriod"];
    settings::generateGrid = settings["generateGrid"];
    settings::from_file = settings["from_file"];
    settings::out_folder = settings["out_folder"];
    settings::unique_folder = settings["uniquefolder"];
    settings::timeLimit = settings["timeLimit"];
}

#include <tuple>



int main() {
  std::cout<<"Running main"<<std::endl;
  try
  {
  loadsettings();
  auto limGrid = std::make_tuple(settings::sizex,settings::sizey , settings::sizez);
  auto limLug = std::make_tuple(settings::minx, settings::maxx, settings::miny, settings::maxy, settings::minz,  settings::maxz);
  auto r = runFlow(limGrid, limLug, settings::generateGrid, settings::from_file ,settings::out_folder,settings::unique_folder);
  r.run(settings::timeLimit, settings::recordPeriod, settings::recordPeriod);
  }
  catch(std::exception& e)
  {
    std::cout<<"Exception happened: "<<e.what()<<std::endl;
  }
  return 0;
}
