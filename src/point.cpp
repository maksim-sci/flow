//
// Created by Danil-Ponkratov
//

#include <point/point.hpp>

#include <tuple>
#include <iostream>

point::point(int x, int y, int z)
  : p(std::make_tuple(x, y, z))
{}

point operator- (const point& p1, const point& p2)
{
  auto [x1, y1, z1] = p1.p;
  auto [x2, y2, z2] = p2.p;
  return point(x1-x2, y1-y2, z1-z2);
}

std::ostream& operator << (std::ostream& os, const point& p)
{
  auto [x, y, z] = p.p;
  os << ' ' << x << ' ' << y << ' ' << z << ' ';
  return os;
}