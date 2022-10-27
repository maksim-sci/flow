//
// Created by Danil-Ponkratov
//

#include <random>
#include <iostream>
#include <ios>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include <types/atom.hpp>
#include <grid/grid.hpp>
#include <params/params.hpp>
#include "types.h"


#include <fmt/format.h>


lug::lug()
{}

lug::lug(int lx, int rx, int ly, int ry, int lz, int rz)
    : lx(lx), rx(rx),
      ly(ly), ry(ry),
      lz(lz), rz(rz)
{}

bool lug::inLug(const atom *a)
{
  auto [x, y, z] = a->p.p;

  return x >= lx && x <= rx && y >= ly && y <= ry && z >= lz && z <= rz;
}

void grid::genElectrode()
{
  for (int x = lx; x < rx; x++)
  {
    for (int y = ly; y < ry; y++)
    {
      atoms[x][y][0] = atom(TypeAtom::ELECTRODE, x, y, 0);
      atoms[x][y][rz-1] = atom(TypeAtom::ELECTRODE, x, y, rz-1);
      cnt_+=2;

    }
  }
}


grid::grid(const int _rx, const int _ry, const int _rz, params* _param, bool generate)
      : rx(_rx)
      , ry(_ry)
      , rz(_rz+2)
      , lx(0)
      , ly(0)
      , lz(0)
      ,param(_param)
{
   //deprecated, will be fixed
   //TODO совместимость с кривоугольными координатами!

  param->E = param->U//С новой версии используется задание напряженности через напряжение между электродами. 
  /distance(pos_t(0,0,0),pos_t(0,0,rz+1))//расстояние между электродами
  *1e+10;//расстояние в ангстремах! поэтому домножаем на размерность, чтобы получить напряженность.
  //хз как быть с кривыми электродами, но что поделать?
  //реально, плак-плак, из-за кривых электродов
  //корректно работает только с прямоугольными координатами, кроме прочего

  std::cout<<"grid initializing"<<std::endl;
  atoms.resize(rx);
  for (auto& y : atoms)
  {
    y.resize(ry);
    for (auto& z : y)
    {
      z.resize(rz);
    }
  }

  for (int x = 0; x < rx; x++)
  {
    for (int y = 0; y < ry; y++)
    {
      for (int z = 0; z < rz; z++)
      {
        atoms[x][y][z] = atom(TypeAtom::EMPTY_INTERSTITIAL_ATOM, x, y, z);
      }
    }
  }
  


  genElectrode();
  
  if (!generate)
  {
    return;
  }
  std::cout<<"Generating grid"<<std::endl;
  static std::random_device rd;
  std::mt19937 gen(rd());
  static std::uniform_int_distribution<> change(0, 10);
  int RIGHT_NUMB = 5;

  for (int x = lx+1; x < rx; x+=2)
  {
    for (int y = ly+1; y < ry; y+=2)
    {
      for (int z = lz+1; z < rz-1; z+=2)
      {
        cnt_++;
        atoms[x][y][z] = atom(TypeAtom::FULL_NODE, x, y, z);
        if (change(gen) == RIGHT_NUMB) {
          atoms[x][y][z] = atom(TypeAtom::VACANCY_NODE_WITH_ELECTRON, x, y, z);
        }
      }
    }
  }

  for (int x = lx; x < rx; x+=2)
  {
    for (int y = ly; y < ry; y+=2)
    {
      for (int z = lz+2; z < rz-1; z+=2)
      {
        if (change(gen) == RIGHT_NUMB) {
          cnt_++;
          atoms[x][y][z] = atom(TypeAtom::INTERSTITIAL_ATOM, x, y, z);
        }
      }
    }
  }

  
}


std::ostream& operator << (std::ostream& os, const grid& g)
{
  for (auto& vx:g.atoms)
  {
    for (auto& vy:vx)
    {
      for ( auto& a:vy)
      {
        if (a.type == TypeAtom::EMPTY_INTERSTITIAL_ATOM)
        {
          continue;
        }
        os << a << std::endl;
      }
    }
  }
  return os;
}

std::istream& operator >> (std::istream& is, grid& g)
{
  g.cnt_ = 0;
  std::string line;
  while (std::getline(is, line))
  {
    if (line.empty())
    {
      continue;
    }

    std::istringstream ss(line);
    atom a;
    ss >> a;
    if (a.type == TypeAtom::ELECTRODE || a.type == TypeAtom::EMPTY_INTERSTITIAL_ATOM)
    {
      continue;
    }

    g.cnt_++;
    auto [ax, ay, az] = a.p.p;
    g.atoms[ax][ay][az] = a;
  }

  return is;
}


void toFile(std::string name, const grid& g)
{
  std::ofstream f(name, std::ios_base::out);
  f << g.cnt_ << std::endl;
  f << getTypesAtom() << std::endl;
  f << g;
}

void fromFile(std::string name, grid& g)
{
  std::ifstream f(name, std::ios_base::in);
  std::string _line;
  std::getline(f, _line); std::getline(f, _line); // строчки с количеством и типами
  f >> g;
}

void grid::calcU()
{
  int lugSize = (rx-lx)*(ry-ly)*(rz-lz);
  for (auto& vx:atoms)
  {
    for(auto&vy:vx)
    {
      for(auto&a:vy)
      {
        {
          auto type = checkElectrodeType(a.p.p);
          bool positiveVoltage = param->U>0;
          if(type==ElectrodeType::POSITIVE)
          {
            if(positiveVoltage)
            {
              a.u = param->U;
            }
            else
            {
              a.u = 0;
            }
            continue;
          }
          else if (type==ElectrodeType::NEGATIVE)
          {
            if(positiveVoltage)
            {
              a.u = 0;
            }
            else
            {
              a.u = param->U;
            }
            continue;
          }
        }
        //capacitor voltage:
        auto [x,y,z] = a.p.p;
        a.u = param->U*(z-lz)/(rz-lz);
        //////////////////
        //"voltage" from lug (pseudo):
        bool lugIsLeft = lug.lz<(rz+lz)/2;
        //TODO support of the right lug;
        int sizeLug = (lug.rz-lug.lz)*(lug.ry-lug.ly)*(lug.rz-lug.lz-1);
        double lug_charge = sizeLug*param->e_charge;
        double distance;
        {
          double dx = dx>=lug.lx?(dx<=lug.rx?0:lug.rx-dx):lug.lx-dx;
          double dy = dy>=lug.ly?(dy<=lug.ry?0:lug.ry-dy):lug.lx-dy;
          double dz = dz>=lug.lz?(dz<=lug.rz?0:lug.rz-dz):lug.lx-dz;
          distance = sqrt(dx*dx+dy*dy+dz*dz)*1e-10;
          //TODO check this!!!
        }
        double lug_effect = (distance<1.1e-6)?0:lug_charge/distance;
        //fmt::print("distance: {} u cond: {}   u effect: {} u result: {}\n",distance, a.u,lug_effect,a.u+lug_effect);
        a.u += lug_effect;
      }
    }
  }
}

atom* grid::get(pos_t pos)
{

  auto [x,y,z] = pos;
  return get(x,y,z);
}

atom* grid::get(int x, int y, int z)
{
  if(lx>x || rx<=x)
  {
    return nullptr;
  }
  if(ly>y || ry<=y)
  {
    return nullptr;
  }
  if(lz>z || rz<=z)
  {
    return nullptr;
  }
  return &atoms[x][y][z];
}

ElectrodeType grid::checkElectrodeType(pos_t pos)
{
  auto* a = get(pos);
  if(a==nullptr) return ElectrodeType::NONE;
  if(a->type!=TypeAtom::ELECTRODE) return ElectrodeType::NONE;
  ElectrodeType electrodeType;
  auto [x,y,z] = pos;
  {
    bool isElectrodeLeft = z<=(rz-lz)/2;
    bool isLeftElectrodePositive = param->U<0;
    bool isElectrodePositive = (- (isElectrodeLeft&&isLeftElectrodePositive)) && (isElectrodeLeft&&isLeftElectrodePositive);
    electrodeType = isElectrodePositive?ElectrodeType::NEGATIVE:ElectrodeType::POSITIVE;
  }
  return electrodeType;
}


void grid::checkElectrode(pos_t pos)
{
  atom* a = get(pos);
  if(a==nullptr) return;
  if(a->type!=TypeAtom::ELECTRODE) return;
  auto [x,y,z] = pos;
  
  ElectrodeType electrodeType = checkElectrodeType(pos);
  for (int dx = 0; dx<=2; dx++)
    for(int dy=0; dy<=2;dy++)
      for(int dz=0; dz<=2;dz++)
      {
        if(dx==0&&dy==0&&dz==0) continue;
        int nx = x+dx,ny = x+dy, nz=x+dz;
        atom* na = get(nx,ny,nz);
        if(na==nullptr) continue;
        na->setElectrodeType(electrodeType);
      }

}


void grid::addLug(int lx_, int rx_, int ly_, int ry_, int lz_, int rz_)
{
  //now lug = additional cubic electrode;
  for (int x=lx_; x<=rx_;x++)
    for (int y=ly_; y<=ry_;y++)
      for (int z=lz_; z<=rz_;z++)
      {

        atoms[x][y][z].setType(TypeAtom::ELECTRODE);
        checkElectrode(pos_t(x,y,z));
      }
  lug = {lx_, rx_, ly_, ry_, lz_, rz_};
}

void grid::updElectrodes()
{
  for(auto& vx:atoms)
    for(auto&vy:vx)
      for(auto& a:vy)
      {
        checkElectrode(a.p.p);
      }
}

void grid::invert()
{
  for(auto& a:atoms)
  {
    for(auto&b:a)
    {
      size_t max = b.size();
      for(size_t cnt = 0; cnt<=max/2;cnt++ )
      {
        std::swap(b[cnt],b[max-cnt]);;
      }
    }
  }
}

  double grid::distance(pos_t first,pos_t second)
  {
    double a = (std::get<0>(second)-std::get<0>(first))*param->a;
    double b = (std::get<1>(second)-std::get<1>(first))*param->b;
    double c = (std::get<2>(second)-std::get<2>(first))*param->c;
    double bb = b+c*param->cosb;
    double cc = c*param->sinb;
    double bbb = bb*param->sina;
    double aa = a+bb*param->cosa;
    return sqrt(aa*aa+bbb*bbb+cc*cc);//TODO remake this?!
  }

 

