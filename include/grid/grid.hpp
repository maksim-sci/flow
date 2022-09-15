//
// Created by Danil-Ponkratov
//

#pragma once

#include <atom/atom.hpp>
#include <electrode/electrode.hpp>
#include "types.h"
#include <params/params.hpp>

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
  grid(const int rx, const int ry, const int rz, params* _param,bool generate);
  void invert();
  friend std::ostream& operator << (std::ostream&, const grid&);
  friend std::istream& operator >> (std::istream&, grid&);

  atom* get(int x, int y, int z);

  void addLug(int lx, int rx, int ly, int ry, int lz, int rz);

  double distance(pos_t first,pos_t second);

  int cnt_ = 0;
  electrode el;
  std::vector< std::vector < std::vector < atom > > > atoms;
  const int rx, ry, rz;
  class lug lug;
  params* param;
};

void toFile(std::string, const grid&);
void fromFile(std::string, grid&);