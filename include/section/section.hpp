//
// Created by Danil-Ponkratov
//

#pragma once

#include <utility>
#include <vector>
#include <ostream>

class section
{
public:
  section();

  void add(const std::pair<int, double>&);
  int get();
  void clear();

  friend std::ostream& operator << (std::ostream&, const section&);

private:
  std::vector<std::pair<int, double>> data;
  double _sum;
};