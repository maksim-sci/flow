//
// Created by Danil-Ponkratov
//

#include <section/section.hpp>

#include <cassert>
#include <random>
#include <utility>
#include <algorithm>
#include <iterator>
#include <ostream>

section::section():_sum(0), data() {}

void section::add(const std::pair<int, double>& p)
{
  assert(p.second <= 1.0);
  data.push_back(p);
  _sum += p.second;
  //printf("seconds %f, curr: %f \n", p.second, _sum);
}

int section::get()
{
  assert(_sum >= 0.99 && _sum <= 1.01);
  assert(data.size() > 0);
  std::random_device rd;
  // перемешаем данные
  auto rngShuffle = std::default_random_engine { rd() };
  std::shuffle(std::begin(data), std::end(data), rngShuffle);

  std::mt19937 rng(rd());
  std::uniform_real_distribution<> dis(0, 1);
  double rn = dis(rng);
  double sum = 0.0;
  auto is = [rn, &sum] (const std::pair<int, double>& el) mutable -> bool { sum = sum + el.second; return sum >= rn; };
  auto res = std::find_if(begin(data), end(data), is);
  return res->first;
}


std::ostream& operator << (std::ostream& os, const section& s)
{
  for (auto& el : s.data)
  {
    os << "{ " << el.first << ' ' << el.second << " } ";
  }
  return os;
}

void section::clear()
{
  data.clear();
  _sum = 0.0;
}
