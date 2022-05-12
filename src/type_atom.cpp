//
// Created by Danil-Ponkratov
//
#include <types/atom.hpp>

#include <iostream>
#include <string>

std::ostream& operator << (std::ostream& os, const TypeAtom& ta)
{
  switch (ta)
  {
  case TypeAtom::FULL_NODE: os << "FULL_NODE"; break;
  case TypeAtom::VACANCY_NODE_WITH_ELECTRON: os << "VACANCY_NODE_WITH_ELECTRON"; break;
  case TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON: os << "VACANCY_NODE_WITHOUT_ELECTRON"; break;
  case TypeAtom::EMPTY_INTERSTITIAL_ATOM: os << "EMPTY_INTERSTITIAL_ATOM"; break;
  case TypeAtom::INTERSTITIAL_ATOM: os << "INTERSTITIAL_ATOM"; break;
  case TypeAtom::ELECTRODE: os << "ELECTRODE"; break;
  }

  return os;
}

std::istream& operator >> (std::istream& is, TypeAtom& ta)
{
  std::string tp;
  is >> tp;
  if (tp == "F" || tp == "FULL_NODE") ta = TypeAtom::FULL_NODE;
  else if (tp == "V+" || tp == "VACANCY_NODE_WITH_ELECTRON") ta = TypeAtom::VACANCY_NODE_WITH_ELECTRON;
  else if (tp == "V0" || tp == "VACANCY_NODE_WITHOUT_ELECTRON") ta = TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON;
  else if (tp == "I0" || tp == "EMPTY_INTERSTITIAL_ATOM") ta = TypeAtom::EMPTY_INTERSTITIAL_ATOM;
  else if (tp == "I" || tp == "INTERSTITIAL_ATOM") ta = TypeAtom::INTERSTITIAL_ATOM;
  else if (tp == "E" || tp == "ELECTRODE") ta = TypeAtom::ELECTRODE;

  return is;
}

std::string shortName(const TypeAtom& ta)
{
  std::string sn;
  switch (ta)
  {
  case TypeAtom::FULL_NODE: sn = "F"; break;
  case TypeAtom::VACANCY_NODE_WITH_ELECTRON: sn = "V+"; break;
  case TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON: sn = "V0"; break;
  case TypeAtom::EMPTY_INTERSTITIAL_ATOM: sn = "I0"; break;
  case TypeAtom::INTERSTITIAL_ATOM: sn = "I"; break;
  case TypeAtom::ELECTRODE: sn = "E"; break;
  }
  return sn;
}

std::string getTypesAtom()
{
  return "F V+ V0 I0 I E";
}