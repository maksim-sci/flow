//
// Created by Danil-Ponkratov
//

#pragma once

#include <types/atom.hpp>
#include <point/point.hpp>

#include <iostream>

class atom
{
public:
  atom();
  atom(TypeAtom typeAtom, int x, int y, int z);

  point p;
  int q;
  TypeAtom type;
};

std::ostream& operator << (std::ostream&, const atom&);
std::istream& operator >> (std::istream&, atom&);

int getQ(TypeAtom type);