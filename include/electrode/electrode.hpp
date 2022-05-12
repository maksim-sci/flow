//
// Created by Danil-Ponkratov
//

#pragma once

#include <atom/atom.hpp>

#include <vector>
#include <array>

class electrode
{
public:
  electrode(int lx, int rx, int ly, int ry, int lz, int rz);

  int getCntAtoms() const;

  friend std::ostream& operator << (std::ostream&, const electrode&);

private:
  const int top = 1;
  const int down = 0;
  std::array< std::vector<atom>, 2> plates;

  int cnt_ = 0;
};

