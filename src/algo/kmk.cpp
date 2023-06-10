#include "kmk.hpp"

#include <algorithm>

namespace algo{

kmk::kmk(grid::Grid *_g)
    : reacts(), recieved_reacts(), calculated_reacts(), sum(0), g(_g){};
std::pair<bool, kmk::react_data> kmk::findAndProcessReact() {
  auto data = chooseReact();
  if (not data.first) {
    recalc();
    data = chooseReact();
  }

  if (not data.first) {
    return data;
  }

  processReact(data.second);

  return data;
};
void kmk::processReact(react_data &data) {
  data.r->Apply(data.f, data.s);
  time+=1/sum;

  auto clear_reacts = [&](auto& atom) mutable{

    //atom-a2 pairs
    {
      auto range = calculated_reacts.equal_range(atom);
      for(auto iter=range.first;iter!=range.second;iter++) {
        auto data = iter->second;
        sum-=data.chance;
      }
      calculated_reacts.erase(range.first,range.second);
    }
    printf("c: %d\n",Count());
    printf("r: %d\n",recieved_reacts.size());
    //a2-atom pairs
    {
      auto range = recieved_reacts.equal_range(atom);
      for(auto iter=range.first;iter!=range.second;iter++) {
        auto& data = iter->second;
        auto range = calculated_reacts.equal_range(data.f);
        auto iter_calculated = range.first;
        while(iter_calculated!=range.second) {
          auto& data = iter_calculated->second;
          if(data.s==atom) {
            sum-=data.chance;
            printf("removed!\n");
            iter_calculated = calculated_reacts.erase(iter_calculated);
          }
          else {
            iter_calculated++;
          }
        }
      }
      recieved_reacts.erase(range.first,range.second);
    }
  };

  clear_reacts(data.f);
  clear_reacts(data.s);


};
std::pair<bool, kmk::react_data> kmk::chooseReact() {

  static std::random_device dev;
  static std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0, sum);

  double rand_targ = dist(rng);
  double top = 0;

  auto iter = calculated_reacts.begin();
  react_data *data = nullptr;
  while (top < sum && iter != calculated_reacts.end()) {
    auto &i_data = iter->second;
    top += i_data.chance;

    data = &i_data;

    iter++;
  }

  if (data != nullptr) {
    return {true, *data};
  }

  return {false, react_data()};
};
void kmk::remove(ptr_react r) { reacts.remove(r); };
void kmk::add(ptr_react r) { reacts.push_back(r); };
void kmk::recalc() {
  recieved_reacts.clear();
  calculated_reacts.clear();
  sum = 0;

  for (auto react : reacts) {
    g->for_each([&](const auto &pos, auto &atom) mutable {
      g->for_each([&](const auto &pos2, auto &atom2) {
        if (react->AreAtomsOk(atom, atom2)) {
          double chance = react->Chance(atom, atom2, (pos2 - pos).abs());

          if (chance > 0) {
              react_data data{atom, atom2, pos, pos2, react, chance};
            calculated_reacts.insert({atom, data});
            recieved_reacts.insert({atom2, data});
            sum += chance;
          }
        }
      });
    });
  }
};

size_t kmk::Count() const { return calculated_reacts.size(); };
double kmk::Sum() const { return sum; };
double kmk::Time() const { return time; };
void kmk::ClearTime() { time = 0; };
} // namespace algo