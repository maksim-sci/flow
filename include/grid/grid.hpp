//
// Created by Danil-Ponkratov
//

#pragma once

#include <atom/atom.hpp>
#include <electrode/electrode.hpp>

#include <vector>
#include <string>

class lug
{
public:
  lug();
  lug(int lx, int rx, int ly, int ry, int lz, int rz);

  bool inLug(const atom* a);

  int lx, rx;
  int ly, ry;
  int lz, rz;
};

class grid
{
public:
  grid(const int rx, const int ry, const int rz, bool generate);

  friend std::ostream& operator << (std::ostream&, const grid&);
  friend std::istream& operator >> (std::istream&, grid&);

  atom* get(int x, int y, int z);

  void addLug(int lx, int rx, int ly, int ry, int lz, int rz);

  int cnt_ = 0;
  electrode el;
  std::vector< std::vector < std::vector < atom > > > atoms;
  const int rx, ry, rz;
  class lug lug;
};

void toFile(std::string, const grid&);
void fromFile(std::string, grid&);