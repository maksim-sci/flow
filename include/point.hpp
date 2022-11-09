//
// Created by Danil-Ponkratov
//

#pragma once

#include "types.h"
#include <tuple>
#include <iostream>

class point
{
public:
  point(int, int, int);

  pos_t p;
};

point operator- (const point&, const point&);
std::ostream& operator << (std::ostream&, const point&);