//
// Created by Danil-Ponkratov
//

#include <types/atom.hpp>
#include <grid/grid.hpp>
#include <electrode/electrode.hpp>

#include <random>
#include <iostream>
#include <ios>
#include <fstream>
#include <string>
#include <sstream>

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

grid::grid(const int rx, const int ry, const int rz, bool generate)
      : el(electrode(0, rx, 0, ry, -1, rz+1))
      , rx(rx)
      , ry(ry)
      , rz(rz)
{
  atoms.resize(rx+1);
  for (auto& y : atoms)
  {
    y.resize(ry+1);
    for (auto& z : y)
    {
      z.resize(rz+1);
    }
  }

  for (int x = 0; x <= rx; x++)
  {
    for (int y = 0; y <= ry; y++)
    {
      for (int z = 0; z <= rz; z++)
      {
        atoms[x][y][z] = atom(TypeAtom::EMPTY_INTERSTITIAL_ATOM, x, y, z);
      }
    }
  }

  if (!generate)
  {
    return;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> change(0, 10);
  int RIGHT_NUMB = 5;

  for (int x = 0; x <= rx; x+=2)
  {
    for (int y = 0; y <= ry; y+=2)
    {
      for (int z = 0; z <= rz; z+=2)
      {
        cnt_++;
        atoms[x][y][z] = atom(TypeAtom::FULL_NODE, x, y, z);
        if (change(gen) == RIGHT_NUMB) {
          atoms[x][y][z] = atom(TypeAtom::VACANCY_NODE_WITH_ELECTRON, x, y, z);
        }
      }
    }
  }

  for (int x = 1; x <= rx; x+=2)
  {
    for (int y = 1; y <= ry; y+=2)
    {
      for (int z = 1; z <= rz; z+=2)
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
  for (int x = 0; x < g.atoms.size(); x++)
  {
    for (int y = 0; y < g.atoms[x].size(); y++)
    {
      for (int z = 0; z < g.atoms[x][y].size(); z++)
      {
        if (g.atoms[x][y][z].type == TypeAtom::EMPTY_INTERSTITIAL_ATOM)
        {
          continue;
        }

        os << g.atoms[x][y][z] << std::endl;
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
    if (a.type == TypeAtom::ELECTRODE)
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
  int cntEl = g.el.getCntAtoms();
  f << g.cnt_ + cntEl << std::endl;
  f << getTypesAtom() << std::endl;
  f << g;
  f << g.el;
}

void fromFile(std::string name, grid& g)
{
  std::ifstream f(name, std::ios_base::in);
  std::string _line;
  std::getline(f, _line); std::getline(f, _line); // строчки с количеством и типами
  f >> g;
}

atom* grid::get(int x, int y, int z)
{
  return &atoms[x][y][z];
}

void grid::addLug(int lx_, int rx_, int ly_, int ry_, int lz_, int rz_)
{
  lug = {lx_, rx_, ly_, ry_, lz_, rz_};
}