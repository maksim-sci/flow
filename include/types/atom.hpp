//
// Created by Danil-Ponkratov
//

#pragma once

#include <iostream>
#include <string>

enum class TypeAtom
{
  FULL_NODE,
  VACANCY_NODE_WITH_ELECTRON,
  VACANCY_NODE_WITHOUT_ELECTRON,
  EMPTY_INTERSTITIAL_ATOM,
  INTERSTITIAL_ATOM,
  ELECTRODE
};

std::ostream& operator << (std::ostream&, const TypeAtom&);
std::istream& operator >> (std::istream&, TypeAtom&);

std::string getTypesAtom();
std::string shortName(const TypeAtom&);