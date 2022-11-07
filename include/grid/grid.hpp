//
// Created by Danil-Ponkratov
//

#pragma once

#include <atom/atom.hpp>
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

  void updElectrodes();
  atom* get(int x, int y, int z);
  atom* get(pos_t pos);
  size_t count() const;
  ElectrodeType checkElectrodeType(pos_t pos);

  void calcU();

  void genElectrode();

  void addLug(int lx, int rx, int ly, int ry, int lz, int rz);

  void checkElectrode(pos_t pos);

  double distance(pos_t first,pos_t second);

  int cnt_ = 0;
  std::vector< std::vector < std::vector < atom > > > atoms;
  const int rx, ry, rz;
  const int lx, ly, lz;
  class lug lug;
  params* param;
};

void toFile(std::string, const grid&);
void fromFile(std::string, grid&);