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
      t(params.E),
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
  E1Shift = {
      point(0,0,-2),point(0,0,2),
      point(0,2,-2),point(0,2,2),
      point(0,-2,-2),point(0,-2,2),
      point(2,-2,-2),point(2,-2,2),
      point(-2,-2,-2),point(-2,-2,2)
  };
  E2Shift = {
      point(0,0,-2),point(0,0,2),
      point(0,2,-2),point(0,2,2),
      point(0,-2,-2),point(0,-2,2),
      point(2,-2,-2),point(2,-2,2),
      point(-2,-2,-2),point(-2,-2,2)
  };
  E3Shift = {
      point(0,0,-2),point(0,0,2),
      point(0,2,-2),point(0,2,2),
      point(0,-2,-2),point(0,-2,2),
      point(2,-2,-2),point(2,-2,2),
      point(-2,-2,-2),point(-2,-2,2)
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
          flow::calc(TypeReaction::E1,a, p->E1mult);
          if(z>=grid->rz-1)
          {
            flow::calcE(TypeReaction::E3,a,p->E3mult);
          }
          break;
        case TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON:
          flow::calc(TypeReaction::R3, a, p->R3mult);
          flow::calc(TypeReaction::R4, a, p->R4mult);
          if(z<=1)
          {
            flow::calcE(TypeReaction::E2, a, p->E2mult);
          }
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
  //std::cout<<sFreq<<" "<< reactionsBox[i].freq<<std::endl;
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
  //std::cout<<"Transition ended"<<std::endl;
}

void flow::transition(reactionData& react)
{

  switch (react.type)
  {
  case TypeReaction::R1:

  //std::cout<<"transitionR1"<<std::endl;

    flow::transitionR1(react);
    statistic[(int) TypeReaction::R1]++;
    break;
  case TypeReaction::R2:
  //std::cout<<"transitionR2"<<std::endl;
    flow::transitionR2(react);
  //std::cout<<"transitionR2 1000-7"<<std::endl;
    statistic[(int) TypeReaction::R2]++;
  //std::cout<<"transitionR2 1000-7"<<std::endl;
    break;
  case TypeReaction::R3:
  //std::cout<<"transitionR3"<<std::endl;
    flow::transitionR3(react);
    statistic[(int) TypeReaction::R3]++;
    break;
  case TypeReaction::R4:
  //std::cout<<"transitionR4"<<std::endl;
    flow::transitionR4(react);
    statistic[(int) TypeReaction::R4]++;
    break;
  case TypeReaction::E1:
  //std::cout<<"transitionE1"<<std::endl;
    flow::transitionE1(react);
    statistic[(int) TypeReaction::E1]++;
    break;
  
  case TypeReaction::E2:
  //std::cout<<"transitionE2"<<std::endl;
    flow::transitionE2(react);
    statistic[(int) TypeReaction::E2]++;
    break;
  
  case TypeReaction::E3:
  //std::cout<<"transitionE3"<<std::endl;
    flow::transitionE3(react);
    statistic[(int) TypeReaction::E3]++;
    break;
  }
  //std::cout<<"889898989"<<std::endl;
}


void flow::calc(TypeReaction rType, class atom* atom, double gain)
{
  std::vector<class point> RShift;
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
  case TypeReaction::E1:
    RShift = E1Shift;
    tp = TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON;
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

    double E = t.getE(atom, sA, gain, grid->lug.inLug(sA));
    //std::cout<<"E: "<<E<<" "<<grid->lug.inLug(sA)<<std::endl;
    double prod = p->y * p->a * p->e * E;
    if (Ea - prod < 0)
    {
      fmt::print("incorrect frequency: reaction type {}, first atom {}, second atom {}\n", rType, *atom, *sA);
      continue;
    }

    double freq = p->v * exp(-( (Ea - prod)/(p->k * p->T) ));
    sFreq += freq;

    reactionData reactionDt= {std::make_tuple(x, y, z), std::make_tuple(xs, ys, zs), freq, rType};
    reactionsBox.push_back(reactionDt);
  }
  //std::cout<<reactionsBox.size()<<std::endl;
}

void flow::transitionR1(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  flow::helperTransition(react.second, TypeAtom::INTERSTITIAL_ATOM);
  grid->cnt_++;

}

void flow::transitionR2(reactionData& react)
{
  //std::cout<<"Transition R2 0mg"<<std::endl;
  flow::helperTransition(react.first, TypeAtom::EMPTY_INTERSTITIAL_ATOM);
  //std::cout<<"Transition R2 0m2"<<std::endl;
  flow::helperTransition(react.second, TypeAtom::INTERSTITIAL_ATOM);
  //std::cout<<"Transition R2 0m3"<<std::endl;
}

void flow::transitionR3(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::FULL_NODE);
  flow::helperTransition(react.second, TypeAtom::EMPTY_INTERSTITIAL_ATOM);
  grid->cnt_--;
}

void flow::transitionR4(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::FULL_NODE);
  flow::helperTransition(react.second, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
}

void flow::helperTransition(std::tuple<int, int, int> c, TypeAtom newType) const
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

void flow::calcE(TypeReaction rType, class atom* atom, double gain)
{
  std::vector<class point> RShift;
  TypeAtom tp;
  double dE;
  double A;
  size_t pos;
  switch (rType) {
  case TypeReaction::E2:
    dE = p->DE2;
    A = p->AE2;
    pos = -1;
    break;
  case TypeReaction::E3:
    dE = p->DE3;
    A = p->AE3;
    pos = grid->rz+1;
    break;
  }

  auto [x,y,z]= atom->p.p;

  auto E = t.getE(atom, grid->get(x,y,pos), gain, grid->lug.inLug(atom));
  double bot = 1-exp(p->e*E*2/(p->kb*p->Temperature));
  if(bot>0) 
  {
    double freq = A*dE * exp(-2/p->le)/exp(p->hconst*(bot));
    if(freq)
    sFreq += freq;
    //std::cout<<freq<<std::endl;
    reactionData reactionDt= {std::make_tuple(x, y, z), std::make_tuple(0, 0, 0), freq, rType};
    reactionsBox.push_back(reactionDt);

  }
}

void flow::transitionE1(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  flow::helperTransition(react.second, TypeAtom::VACANCY_NODE_WITH_ELECTRON);
}
void flow::transitionE2(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
}
void flow::transitionE3(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITH_ELECTRON);
  I++;

  //std::cout<<"E3"<<std::endl;
}

reactionData::reactionData(std::tuple<int, int, int>&& first, std::tuple<int, int, int>&& second, double freq, TypeReaction type)
    :   first(first),
      second(second),
      freq(freq),
      type(type)
{}


bool filter(int x, int x_, int lm)
{
  return x + x_ >= 0 && x + x_ <= lm;
}