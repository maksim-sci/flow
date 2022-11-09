//
// Created by Danil-Ponkratov
//

#pragma once

#include <atom.hpp>
#include <grid.hpp>

class tension
{
public:
  explicit tension(double E);

  double getE(const atom* first, const atom* second, double gain, bool isLug);

private:
  double E;
};
