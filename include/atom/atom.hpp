//
// Created by Danil-Ponkratov
//

#pragma once

#include <types/atom.hpp>
#include <point/point.hpp>

#include <iostream>
enum class ElectrodeType
{
  NONE,
  POSITIVE,
  NEGATIVE
};

class atom
{
public:
  atom();
  atom(TypeAtom typeAtom, int x, int y, int z);
  void setElectrodeType(ElectrodeType type);
  void setType(TypeAtom type);
  bool isVacancy();

  point p;
  int q;
  double u;
  ElectrodeType meta;
  TypeAtom type;
};

std::ostream& operator << (std::ostream&, const atom&);
std::istream& operator >> (std::istream&, atom&);

int getQ(TypeAtom type);