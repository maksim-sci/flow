#include "kmk.hpp"
#include "geometry/vector.hpp"
#include "grid/atom/atom.hpp"

#include <algorithm>

namespace algo {

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
  if (apply) {
    data.r->Apply(data.f, data.s);
  }
  time += 1 / sum;

  auto clear_reacts = [&](auto &atom) mutable {
    // atom-a2 pairs
    {
      auto range = calculated_reacts.equal_range(atom);
      for (auto iter = range.first; iter != range.second; iter++) {
        auto data = iter->second;
        sum -= data.chance;
      }
      calculated_reacts.erase(range.first, range.second);
    }
    // a2-atom pairs
    {
      auto range = recieved_reacts.equal_range(atom);
      for (auto iter = range.first; iter != range.second; iter++) {
        auto &data = iter->second;
        auto range = calculated_reacts.equal_range(data.f);
        auto iter_calculated = range.first;
        while (iter_calculated != range.second) {
          auto &data = iter_calculated->second;
          if (data.s == atom) {
            sum -= data.chance;
            iter_calculated = calculated_reacts.erase(iter_calculated);
          } else {
            iter_calculated++;
          }
        }
      }
      recieved_reacts.erase(range.first, range.second);
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
  while (top < rand_targ && iter != calculated_reacts.end()) {
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
void kmk::remove(ptr_react r) {

  reacts.erase(std::remove_if(reacts.begin(), reacts.end(),
                              [&r](auto &r2) { return r2 == r; }));
};
void kmk::add(ptr_react r) { reacts.push_back(r); };
void kmk::recalc() {
  recieved_reacts.clear();
  calculated_reacts.clear();
  sum = 0;

  if(cache) {
    for(auto&cachedData:cacheData) {
      if (cachedData.r->AreAtomsOk(cachedData.f, cachedData.s)) {
        double chance = cachedData.r->Chance(cachedData.f, cachedData.s, (cachedData.fp - cachedData.sp).abs());

        if (chance > 0) {
          react_data data{cachedData.f, cachedData.s, cachedData.fp, cachedData.sp, cachedData.r, chance};
          calculated_reacts.insert({data.f, data});
          recieved_reacts.insert({data.s, data});
          sum += chance;
        }
      }
    }
    return;
  }

  for (auto react : reacts) {

    geometry::Vector pos1;
    std::shared_ptr<grid::atom::Atom> atom1;

    auto atom_enum = [&pos1, &atom1, &react, this](const auto &pos2,
                                                   auto &atom2) {
      if (react->AreAtomsOk(atom1, atom2)) {
        double chance = react->Chance(atom1, atom2, (pos2 - pos1).abs());

        if (chance > 0) {
          react_data data{atom1, atom2, pos1, pos2, react, chance};
          calculated_reacts.insert({atom1, data});
          recieved_reacts.insert({atom2, data});
          sum += chance;
        }
      }
    };

    g->for_each([&atom_enum, &pos1, &atom1, &react, this](const auto &pos,
                                                          auto &atom) mutable {
      pos1 = pos;
      atom1 = atom;
      g->for_each(pos1, react->Distance(), atom_enum);
    });
  }
};

size_t kmk::Count() const { return calculated_reacts.size(); };
double kmk::Sum() const { return sum; };
double kmk::Time() const { return time; };
void kmk::ClearTime() { time = 0; };

void kmk::Cache(bool mode) {
  cache = mode;

  cacheData.clear();

  if (mode) {
    for (auto react : reacts) {

      geometry::Vector pos1;
      std::shared_ptr<grid::atom::Atom> atom1;

      auto atom_enum = [&pos1, &atom1, &react, this](const auto &pos2,
                                                     auto &atom2) {
        double distance = (pos1-pos2).abs();
        if (react->Distance()>distance) {
          cacheData.push_back({atom1,atom2,react,pos1,pos2});
        }
      };

      g->for_each([&atom_enum, &pos1, &atom1, &react,
                   this](const auto &pos, auto &atom) mutable {
        pos1 = pos;
        atom1 = atom;
        g->for_each(pos1, react->Distance(), atom_enum);
      });
    }
  }
}
} // namespace algo