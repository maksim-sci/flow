//
// Created by Danil-Ponkratov
//

#include <atom/atom.hpp>
#include <types/atom.hpp>
#include <electrode/electrode.hpp>

#include <vector>
#include <array>
#include <iostream>

electrode::electrode(int lx, int rx, int ly, int ry, int lz, int rz)
{
  plates = std::array< std::vector<atom>, 2 >{};
  for (int x = lx; x <= rx; x++)
  {
    for (int y = ly; y <= ry; y++)
    {
      plates[down].push_back(atom(TypeAtom::ELECTRODE, x, y, lz));
      plates[top].push_back(atom(TypeAtom::ELECTRODE, x, y, rz));
      cnt_+=2;
    }
  }
}

int electrode::getCntAtoms() const
{
  return cnt_;
}

std::ostream& operator << (std::ostream& os, const electrode& e)
{
  int sz = e.plates[e.down].size();
  for (int i = 0; i < sz; i++)
  {
    os << e.plates[e.down][i] << std::endl;
    os << e.plates[e.top][i] << std::endl;
  }
  return os;
}

