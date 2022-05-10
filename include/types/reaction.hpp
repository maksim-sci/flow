//
// Created by Danil-Ponkratov
//

#pragma once

#include <iostream>

enum class TypeReaction
{
  R1,
  R2,
  R3,
  R4,
  E1,
  E2,
  E3
};

std::ostream& operator << (std::ostream&, const TypeReaction&);
