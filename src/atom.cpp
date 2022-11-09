//
// Created by Danil-Ponkratov
//

#include <types/atom.hpp>
#include <point.hpp>
#include <atom.hpp>

#include <iostream>

atom::atom(TypeAtom typeAtom, int x, int y, int z)
    : p(point(x, y, z)),
      type(typeAtom),
      meta(ElectrodeType::NONE),
      u(0)
{
  q = getQ(typeAtom);
}


atom::atom()
    : p(point(-100500, -100500, -100500)),
      type(TypeAtom::EMPTY_INTERSTITIAL_ATOM),
      q(0),
      u(0)
{}

std::ostream& operator << (std::ostream& os, const atom& a)
{
  os << shortName(a.type);
  os << a.p;
  return os;
}

std::istream& operator >> (std::istream& is, atom& a)
{
  int x; int y; int z; TypeAtom ta;
  is >> ta >> x >> y >> z;
  a = atom(ta, x, y, z);
  return is;
}

bool atom::isVacancy()
{
  return (type==TypeAtom::VACANCY_NODE_WITH_ELECTRON || type==TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
}

void atom::setType(TypeAtom _type)
{
  if(type==_type) return;
  bool wasVacancy = isVacancy();
  type = _type;
  if((!isVacancy())&&wasVacancy)
  {
    meta = ElectrodeType::NONE;
  }
}


void atom::setElectrodeType(ElectrodeType _type)
{
  meta = _type;
  if(!isVacancy()) return;
  if(type==TypeAtom::VACANCY_NODE_WITH_ELECTRON && _type==ElectrodeType::POSITIVE)
  {
    setType(TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
    q = getQ(type);
  }
  else if(type==TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON && _type==ElectrodeType::NEGATIVE)
  {
    setType(TypeAtom::VACANCY_NODE_WITH_ELECTRON);
    q = getQ(type);
  }
}


int getQ(TypeAtom type)
{
  int q;
  switch (type)
  {
  case TypeAtom::FULL_NODE: q = 0; break;
  case TypeAtom::VACANCY_NODE_WITH_ELECTRON: q = 1; break;
  case TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON: q = 0; break;
  case TypeAtom::EMPTY_INTERSTITIAL_ATOM: q = 0; break;
  case TypeAtom::INTERSTITIAL_ATOM: q = -1; break;
  case TypeAtom::ELECTRODE: q = 0; break;
  }

  return q;
}