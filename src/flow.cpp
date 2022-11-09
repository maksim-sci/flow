//
// Created by Danil-Ponkratov
//

#include <tuple>
#include <cmath>
#include <vector>
#include <iostream>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <flow.hpp>
#include <grid.hpp>
#include <params.hpp>
#include <section.hpp>
#include <types/atom.hpp>
#include <point.hpp>
#include <atom.hpp>

#include <sgs.hpp>

constexpr double min_valuable_freq = 1e-10;

flow::flow(class grid &grid, class section &section, class params &params)
    : grid(&grid),
      section(&section),
      p(&params),
      tension(params.E),
      dq(0),
      dt(INFINITY),
      statistic(),
      bcalcu(true)
{
  std::cout << "flow initializing" << std::endl;
  R1Shift = {
      pos_t(0, 1, 1), pos_t(0, -1, 1),
      pos_t(0, 1, -1), pos_t(0, -1, -1),
      pos_t(1, 1, 0), pos_t(1, -1, 0),
      pos_t(-1, 1, 0), pos_t(-1, -1, 0)};
  R2Shift = {
      pos_t(-2, 0, 0), pos_t(0, 2, 0),
      pos_t(2, 0, 0), pos_t(0, -2, 0),
      pos_t(0, 0, 2), pos_t(0, 0, -2)};
  R3Shift = {
      pos_t(-1, -1, -1), pos_t(-1, 1, -1),
      pos_t(1, -1, -1), pos_t(1, 1, -1),
      pos_t(-1, -1, 1), pos_t(-1, 1, 1),
      pos_t(1, -1, 1), pos_t(1, 1, 1)};
  R4Shift = {
      pos_t(-2, 0, 0), pos_t(0, 2, 0),
      pos_t(2, 0, 0), pos_t(0, -2, 0),
      pos_t(0, 0, 2), pos_t(0, 0, -2)};

  E1Shift = {
      pos_t(-2, 0, 0), pos_t(0, 2, 0),
      pos_t(2, 0, 0), pos_t(0, -2, 0),
      pos_t(0, 0, 2), pos_t(0, 0, -2)};

  E2Shift = {
      pos_t(0, 0, 1), pos_t(0, 0, -1)};

  E3Shift = {
      pos_t(0, 0, 1), pos_t(0, 0, -1)};
  statistic = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  clearStats();
  reactionsBox.reserve(1024);
}

void flow::clearStats()
{
  statistic[(int)TypeReaction::R1] = 0;
  statistic[(int)TypeReaction::R2] = 0;
  statistic[(int)TypeReaction::R3] = 0;
  statistic[(int)TypeReaction::R4] = 0;
  statistic[(int)TypeReaction::E1] = 0;
  statistic[(int)TypeReaction::E2] = 0;
  statistic[(int)TypeReaction::E3] = 0;
}

void flow::step()
{
  if(bcalcu)
  {
    grid->calcU();
    bcalcu=false;
  }
  reactionsBox.clear();
  sFreq = 0;
  for (auto& vx:grid->atoms)
  {
    for (auto& vy:vx)
    {
      for (auto& a:vy)
      {
        switch (a.type)
        {
        case TypeAtom::FULL_NODE:
          flow::calc(TypeReaction::R1, &a, p->R1mult);
          break;
        case TypeAtom::VACANCY_NODE_WITH_ELECTRON:
          flow::calc(TypeReaction::E1,&a, p->E1mult);
          //flow::calcE(TypeReaction::E3,a,p->E3mult);
          break;
        case TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON:
          flow::calc(TypeReaction::R3, &a, p->R3mult);
          flow::calc(TypeReaction::R4, &a, p->R4mult);
          //flow::calcE(TypeReaction::E2, a, p->E2mult);
          break;
        case TypeAtom::EMPTY_INTERSTITIAL_ATOM:
          break;
        case TypeAtom::INTERSTITIAL_ATOM:
          flow::calc(TypeReaction::R2, &a, p->R2mult);
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
    fmt::print("{} {}: happened assert in section->get()\n", __FILE__, __LINE__);
    return;
  }

  reactionData react = reactionsBox[idx];
  flow::transition(react);
}

void flow::transition(reactionData &react)
{
  dt += 1 / react.freq;
  switch (react.type)
  {
  case TypeReaction::R1:

    flow::transitionR1(react);
    statistic[(int)TypeReaction::R1]++;
    bcalcu = true;
    break;
  case TypeReaction::R2:

    flow::transitionR2(react);
    statistic[(int)TypeReaction::R2]++;
    bcalcu = true;
    break;
  case TypeReaction::R3:

    flow::transitionR3(react);
    statistic[(int)TypeReaction::R3]++;
    bcalcu = true;
    break;
  case TypeReaction::R4:

    flow::transitionR4(react);
    statistic[(int)TypeReaction::R4]++;
    bcalcu = true;
    break;
  case TypeReaction::E1:
    flow::transitionE1(react);
    statistic[(int)TypeReaction::E1]++;
    bcalcu = true;
    break;

  case TypeReaction::E2:
    flow::transitionE2(react);
    statistic[(int)TypeReaction::E2]++;
    bcalcu = true;
    break;

  case TypeReaction::E3:
    flow::transitionE3(react);
    statistic[(int)TypeReaction::E3]++;
    bcalcu = true;
    break;
  }
  checkElectrode(react.first, react.freq);
  checkElectrode(react.second, react.freq);
}

void flow::calc(TypeReaction rType, class atom *atom, double gain)
{
  std::vector<pos_t> RShift;
  TypeAtom tp;
  double Ea;
  TypeAtom check;
  switch (rType)
  {
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
    Ea = p->ER1;
    break;
  case TypeReaction::E2:
    RShift = E2Shift;
    tp = TypeAtom::ELECTRODE;
    Ea = p->ER1;
    break;
  case TypeReaction::E3:
    RShift = E3Shift;
    tp = TypeAtom::ELECTRODE;
    Ea = p->ER1;
    break;
  }

  for (auto &sh : RShift)
  {
    auto [x_, y_, z_] = sh;
    auto [x, y, z] = atom->p.p;
    pos_t pos2 = pos_t(x + x_, y + y_, z + z_);
    auto [x2,y2,z2] = pos2;
    class atom *sA = grid->get(pos2);

    if (sA==nullptr)
      continue;

    auto a = grid->get(x_,y_,z_);

    if(a==nullptr)
     continue;

    if (sA->type != tp)
      continue;
    double prod = p->e * (sA->u-a->u); //1e+10 - distance must be in angstrom //TODO - fix this!!!
    // if (Ea - prod < 0)
    // {
    //   //idk why there was this msg
    //   //fmt::print("{} {}: incorrect frequency: reaction type {}, first atom {}, second atom {}\n", __FILE__, __LINE__, rType, *atom, *sA);
    //   continue;
    // }
    //?!

    double freq = p->v * exp(-((Ea) / (BOLZMAN * p->T)));
    //fmt::print("E: {} Ea: {} prod: {} freq: {} sA->u-a->u: {}\n",a->u*p->E,Ea,prod,freq,sA->u-a->u);
    if(freq<=min_valuable_freq)
    {
      continue;
    }
    sFreq += freq;
    
    reactionData reactionDt = {std::make_tuple(x, y, z), std::make_tuple(x2,y2,z2), freq, rType};
    reactionsBox.push_back(reactionDt);
  }
}

void flow::calcE(TypeReaction rType, class atom* atom, double gain)
{
  std::vector<pos_t> RShift;
  TypeAtom tp;
  double dE;
  double A;
  size_t pos;
  switch (rType) {
  case TypeReaction::E2:
    dE = p->DE2;
    A = p->AE2;
    RShift = E2Shift;
    break;
  case TypeReaction::E3:
    dE = p->DE3;
    A = p->AE3;
    RShift = E3Shift;
    break;
  }

  auto [x,y,z]= atom->p.p;
  for (auto& pp:RShift)
  {
    auto[dx,dy,dz]= pp;
    auto [x2,y2,z2] = pos_t(x+dx,y+dy,z+dz);
    auto* a2 =grid->get(x2,y2,z2);
    if(a2->type==TypeAtom::ELECTRODE)
    {
      auto E = tension.getE(atom, a2, gain, grid->lug.inLug(atom));
      double bot = 1-exp(p->e*E*2/(BOLZMAN*p->Temperature));
      if(bot>0) 
      {
        double freq = A*dE * exp(-2/p->le)/exp(PLANCK*(bot));
        if(freq<=min_valuable_freq)
        {
          continue;
        }
        sFreq += freq;
        reactionData reactionDt= {std::make_tuple(x, y, z), std::make_tuple(x2, y2, z2), freq, rType};
        reactionsBox.push_back(reactionDt);

      }
    }
  }
}

void flow::transitionR1(reactionData &react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  flow::helperTransition(react.second, TypeAtom::INTERSTITIAL_ATOM);
  grid->cnt_++;
}

void flow::transitionR2(reactionData &react)
{
  flow::helperTransition(react.first, TypeAtom::EMPTY_INTERSTITIAL_ATOM);
  flow::helperTransition(react.second, TypeAtom::INTERSTITIAL_ATOM);
}

void flow::transitionR3(reactionData &react)
{
  flow::helperTransition(react.first, TypeAtom::FULL_NODE);
  flow::helperTransition(react.second, TypeAtom::EMPTY_INTERSTITIAL_ATOM);
  grid->cnt_--;
}

void flow::transitionR4(reactionData &react)
{
  flow::helperTransition(react.first, TypeAtom::FULL_NODE);
  flow::helperTransition(react.second, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
}

void flow::transitionE1(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  flow::helperTransition(react.second, TypeAtom::VACANCY_NODE_WITH_ELECTRON);
}
void flow::transitionE2(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITH_ELECTRON);
}
void flow::transitionE3(reactionData& react)
{
  flow::helperTransition(react.first, TypeAtom::VACANCY_NODE_WITHOUT_ELECTRON);
  //dq++;

  //std::cout<<"E3"<<std::endl;
}

void flow::checkElectrode(pos_t c, double freq)
{
  atom* a = grid->get(c);
  if(a==nullptr) return;
  grid->checkElectrode(c);
}

void flow::helperTransition(pos_t c, TypeAtom newType) const
{
  auto [x, y, z] = c;
  atom *a = grid->get(x, y, z);
  a->setType(newType);
  a->q = getQ(newType);
}

void flow::getStatistic()
{
  std::cout << "R1: " << statistic[(int)TypeReaction::R1]
            << " R2: " << statistic[(int)TypeReaction::R2]
            << " R3: " << statistic[(int)TypeReaction::R3]
            << " R4: " << statistic[(int)TypeReaction::R4]
            << " E1: " << statistic[(int)TypeReaction::E1]
            << " E2: " << statistic[(int)TypeReaction::E2]
            << " E3: " << statistic[(int)TypeReaction::E3]
            << std::endl;
}

reactionData::reactionData(pos_t &&first, pos_t &&second, double freq, TypeReaction type)
    : first(first),
      second(second),
      freq(freq),
      type(type)
{
}

bool filter(int x, int x_, int lm)
{
  return x + x_ >= 0 && x + x_ <= lm;
}

double flow::getI()
{
  return dq*ELCHARGE/dt;
}
void flow::clearI()
{
  dq = 0;
  dt = 0;
}