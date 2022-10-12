//
// Created by Danil-Ponkratov
//

#pragma once

class params
{
public:
  params();
  explicit params(double E);

  const double Ea1, Ea2, Ea3, Ea4, ER1, ER2, ER3, y, T, k, e, v;
  double E;
  double R1mult;
  double R2mult;
  double R3mult;
  double R4mult;
  double E1mult;
  double E2mult;
  double E3mult;
  double AE1;
  double AE2;
  double AE3;
  double DE2;
  double DE3;
  double le;
  double hconst;
  double kb;
  double Temperature;
  double a;
  double b;
  double c;
  double alpha;
  double beta;
  double gamma;
  double sina;
  double sinb;
  double sinc;
  double cosa;
  double cosb;
  double cosc;
  double U;
  double e_charge;
  double vac_size;
};