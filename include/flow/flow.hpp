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
  void transition(reactionData& react);

  void transitionR1(reactionData& react);
  void transitionR2(reactionData& react);
  void transitionR3(reactionData& react);
  void transitionR4(reactionData& react);
  void helperTransition(pos_t c, TypeAtom newType) const;

  void getStatistic();

  void calcE(TypeReaction rType, class atom* atom, double gain);
  grid* grid;
  section* section;
  params* p;
  tension t;
  double sFreq = 0.0;

  std::vector<class point> R1Shift;
  std::vector<class point> R2Shift;
  std::vector<class point> R3Shift;
  std::vector<class point> R4Shift;
  std::vector<class point> E1Shift;
  std::vector<class point> E2Shift;
  std::vector<class point> E3Shift;
  std::vector<int> statistic;
  std::vector<reactionData> reactionsBox;

  double I;
};

bool filter(int x, int x_, int lm);