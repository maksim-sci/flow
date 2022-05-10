//
// Created by Danil-Ponkratov
//

#include <types/reaction.hpp>

#include <iostream>

std::ostream& operator << (std::ostream& os, const TypeReaction& r)
{
  switch (r)
  {
  case TypeReaction::R1: os << "R1"; break;
  case TypeReaction::R2: os << "R2"; break;
  case TypeReaction::R3: os << "R3"; break;
  case TypeReaction::R4: os << "R4"; break;
  case TypeReaction::E1: os << "E1"; break;
  case TypeReaction::E2: os << "E2"; break;
  case TypeReaction::E3: os << "E3"; break;
  }

  return os;
}