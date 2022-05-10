//
// Created by Danil-Ponkratov
//

#pragma once

#include <atom/atom.hpp>
#include <grid/grid.hpp>

class tension
{
public:
  explicit tension(double E);

  double getE(const atom* first, const atom* second, double gain, bool isLug);

private:
  double E;
};
