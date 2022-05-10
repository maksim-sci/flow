//
// Created by Danil-Ponkratov
//

#pragma once

#include <tuple>
#include <iostream>

class point
{
public:
  point(int, int, int);

  std::tuple<int, int, int> p;
};

point operator- (const point&, const point&);
std::ostream& operator << (std::ostream&, const point&);