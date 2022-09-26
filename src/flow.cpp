//
// Created by Danil-Ponkratov
//

#include <tuple>
#include <cmath>
#include <vector>
#include <iostream>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <flow/flow.hpp>
#include <grid/grid.hpp>
#include <params/params.hpp>
#include <section/section.hpp>
#include <types/atom.hpp>
#include <point/point.hpp>




flow::flow(class grid& grid, class section& section, class params& params)
  :   grid(&grid),
      section(&section),
      p(&params),
      tension(params.E),
      I(0)
{
  R1Shift = {
      point(0, 1, 1),   point(0, -1, 1),
      point(0, 1, -1),  point(0, -1, -1),
      point(1, 1, 0),   point(1, -1, 0),
      point(-1, 1, 0),  point(-1, -1, 0)
  };
  R2Shift = {
      point(-2, 0, 0),  point(0, 2, 0),
      point(2, 0, 0),   point(0, -2, 0),
      point(0, 0, 2),   point(0, 0, -2)
  };
  R3Shift = {
      point(-1, -1, -1), point(-1, 1, -1),
      point(1, -1, -1),  point(1, 1, -1),
      point(-1, -1, 1),  point(-1, 1, 1),
      point(1, -1, 1),   point(1, 1, 1)
  };
  R4Shift = {
      point(-2, 0, 0),  point(0, 2, 0),
      point(2, 0, 0),   point(0, -2, 0),
      point(0, 0, 2),   point(0, 0, -2)
  };
  statistic = {0,0,0,0,0,0,0,0,0,0,0};
  reactionsBox.reserve(1024);
}

void flow::run()
{
  reactionsBox.clear();
  sFreq = 0;
  for (int x = 0; x < grid->atoms.size(); x++)
  {
    for (int y = 0; y < grid->atoms[x].size(); y++)
    {
      for (int z = 0; z < grid->atoms[x][y].size(); z++)
      {
        atom* a = grid->get(x, y, z);
        switch (a->type)
        {
        case TypeAtom::FULL_NODE:
          flow::calc(TypeReaction::R1, a, p->R1mult);
          break;
        case TypeAtom::VACANCY_NODE_WITH_ELECTRON:
          break;
        case TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON:
          flow::calc(TypeReaction::R3, a, p->R3mult);
          flow::calc(TypeReaction::R4, a, p->R4mult);
          break;
        case TypeAtom::EMPTY_INTERSTITIAL_ATOM:
          break;
        case TypeAtom::INTERSTITIAL_ATOM:
          flow::calc(TypeReaction::R2, a, p->R2mult);
          break;
        case TypeAtom::ELECTRODE:
          break;
        }
      }
    }
  }
  int i;
  for (i = 0; i < reactionsBox.size(); i++)
  {
    double freq = reactionsBox[i].freq / sFreq;
    section->add({i, freq});
  }
  int idx = section->get();
  if (idx < 0 || idx > reactionsBox.size())
  {
    fmt::print("happened assert in section->get()\n");
    return;
  }
  
  reactionData react = reactionsBox[idx];
  flow::transition(react);
}

void flow::transition(reactionData& react)
{

  switch (react.type)
  {
  case TypeReaction::R1:


    flow::transitionR1(react);
    statistic[(int) TypeReaction::R1]++;
    break;
  case TypeReaction::R2:
  
    flow::transitionR2(react);
  
    statistic[(int) TypeReaction::R2]++;
  
    break;
  case TypeReaction::R3:
  
    flow::transitionR3(react);
    statistic[(int) TypeReaction::R3]++;
    break;
  case TypeReaction::R4:
  
    flow::transitionR4(react);
    statistic[(int) TypeReaction::R4]++;
    break;
  }
  
}


void flow::calc(TypeReaction rType, class atom* atom, double gain)
{
  static std::vector<class point> RShift;
  TypeAtom tp;
  double Ea;

  switch (rType) {
  case TypeReaction::R1:
    RShift = R1Shift;
    tp = TypeAtom::EMPTY_INTERSTITIAL_ATOM;
    Ea = p->Ea1;
    break;
  case TypeReaction::R2:
    RShift = R2Shift;
    tp = TypeAtom::EMPTY_INTERSTITIAL_ATOM;
    Ea = p->Ea2;
    break;
  case TypeReaction::R3:
    RShift = R3Shift;
    tp = TypeAtom::INTERSTITIAL_ATOM;
    Ea = p->Ea3;
    break;
  case TypeReaction::R4:
    RShift = R4Shift;
    tp = TypeAtom::FULL_NODE;
    Ea = p->Ea4;
    break;
  }

  for (auto& sh : RShift)
  {
    auto [x_, y_, z_] = sh.p;
    auto [x, y, z] = atom->p.p;
    if ( (! filter(x, x_, grid->rx)) || (! filter(y, y_, grid->ry)) ||
        (! filter(z, z_, grid->rz)) )
      continue;

    auto [xs, ys, zs] = std::make_tuple(x + x_, y + y_, z + z_);
    class atom* sA = grid->get(xs, ys, zs);
    if (sA->type != tp)
      continue;

    double E = tension.getE(atom, sA, gain, grid->lug.inLug(sA));
    
    double prod = p->y * p->a * p->e * E;
    if (Ea - prod < 0)
    {
      fmt::print("incorrect frequency: reaction type {}, first atom {}, second atom {}\n", rType, *atom, *sA);
      continue;
    }

    double freq = p->v * exp(-( (Ea - prod)/(p->k * p->T) ));
    sFreq += freq;
    // if(rType==TypeReaction::R3)
    // {
    //   std::cout<<freq<<std::endl;
    // }
    reactionData reactionDt= {std::make_tuple(x, y, z), std::make_tuple(xs, ys, zs), freq, rType};
    reactionsBox.push_back(reactionDt);
  }
  
}

void flow::transitionR1(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  flow::helperTransition(react.second, TypeAtom::INTERSTITIAL_ATOM);
  I+=p->e*2*grid->distance(react.first,react.second)*react.freq;
  grid->cnt_++;

}

void flow::transitionR2(reactionData& react)
{
  
  flow::helperTransition(react.first, TypeAtom::EMPTY_INTERSTITIAL_ATOM);
  flow::helperTransition(react.second, TypeAtom::INTERSTITIAL_ATOM);
  I+=p->e*react.freq;

  
}

void flow::transitionR3(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::FULL_NODE);
  flow::helperTransition(react.second, TypeAtom::EMPTY_INTERSTITIAL_ATOM);
  grid->cnt_--;
  I+=p->e*(-2)*react.freq;

}

void flow::transitionR4(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::FULL_NODE);
  flow::helperTransition(react.second, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  I+=p->e*(-1)*react.freq;

}

void flow::helperTransition(pos_t c, TypeAtom newType) const
{
  auto [x, y, z] = c;
  atom* a = grid->get(x, y, z);
  a->type = newType;
  a->q = getQ(newType);
}

void flow::getStatistic()
{
  std::cout << "R1: "  << statistic[(int) TypeReaction::R1]
            << " R2: " << statistic[(int) TypeReaction::R2]
            << " R3: " << statistic[(int) TypeReaction::R3]
            << " R4: " << statistic[(int) TypeReaction::R4]
            << std::endl;
}


reactionData::reactionData(pos_t&& first, pos_t&& second, double freq, TypeReaction type)
    :   first(first),
      second(second),
      freq(freq),
      type(type)
{}


bool filter(int x, int x_, int lm)
{
  return x + x_ >= 0 && x + x_ <= lm;
}