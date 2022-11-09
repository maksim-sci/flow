//
// Created by Danil-Ponkratov
//

#include <tension.hpp>

tension::tension(double E)
    : E(E)
{}

double tension::getE(const atom* first, const atom* second, double gain, bool isLug)
{
  double ret = E;
  if (isLug)
    ret = gain * ret;

  point sub = second->p - first->p;
  auto [x, y, z] = sub.p;
  if (z < 0)
    ret = -ret;

  return ret;
}
