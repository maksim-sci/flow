//
// Created by Danil-Ponkratov
//

#pragma once

#include "types.h"
#include <grid/grid.hpp>
#include <section/section.hpp>
#include <params/params.hpp>
#include <atom/atom.hpp>
#include <point/point.hpp>
#include <tension/tension.hpp>

#include <vector>
#include <types/reaction.hpp>

class reactionData
{
public:
  reactionData(pos_t&& first, pos_t&& second, double freq, TypeReaction type);

  pos_t first;
  pos_t second;
  double freq;
  TypeReaction type;
};

class flow
{
public:
  flow(class grid& grid, class section& section, class params& params);

  void run();
  void calc(TypeReaction rType, class atom* atom, double gain);
  void calcE(TypeReaction rType, class atom* atom, double gain);
  void calcE1(TypeReaction rType, class atom* atom, double gain);
  void transition(reactionData& react);

  void transitionR1(reactionData& react);
  void transitionR2(reactionData& react);
  void transitionR3(reactionData& react);
  void transitionR4(reactionData& react);
  void transitionE1(reactionData& react);
  void transitionE2(reactionData& react);
  void transitionE3(reactionData& react);

  void helperTransition(pos_t c, TypeAtom newType) const;
  void checkElectrode(pos_t c, double freq);

  void clearStats();
  void getStatistic();

  
  grid* grid;
  section* section;
  params* p;
  tension tension;
  double sFreq = 0.0;

  std::vector<pos_t> R1Shift;
  std::vector<pos_t> R2Shift;
  std::vector<pos_t> R3Shift;
  std::vector<pos_t> R4Shift;
  std::vector<pos_t> E1Shift;
  std::vector<pos_t> E2Shift;
  std::vector<pos_t> E3Shift;
  
  std::vector<int> statistic;
  std::vector<reactionData> reactionsBox;

  size_t dq;
  double dt;

  double getI();
  void clearI();
};

bool filter(int x, int x_, int lm);