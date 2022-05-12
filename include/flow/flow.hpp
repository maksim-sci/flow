//
// Created by Danil-Ponkratov
//

#pragma once

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
  reactionData(std::tuple<int, int, int>&& first, std::tuple<int, int, int>&& second, double freq, TypeReaction type);

  std::tuple<int, int, int> first;
  std::tuple<int, int, int> second;
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
  void helperTransition(std::tuple<int, int, int> c, TypeAtom newType) const;

  void getStatistic();

  grid* grid;
  section* section;
  params* p;
  tension t;
  double sFreq = 0.0;

  std::vector<class point> R1Shift;
  std::vector<class point> R2Shift;
  std::vector<class point> R3Shift;
  std::vector<class point> R4Shift;
  std::array<int, 4> statistic;
  std::vector<reactionData> reactionsBox;
};

bool filter(int x, int x_, int lm);