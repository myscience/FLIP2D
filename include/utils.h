#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <random>
#include <assert.h>

#define INVALID UINT32_MAX
#define KNOWN 0
#define UNKWN (UINT32_MAX - 1)

#define FLUID 0
#define EMPTY 1
#define SOLID 2


class Grid;
class VectorGrid;

struct Rand {
  std::random_device rd;
  std::ranlux48_base rg;

  std::uniform_real_distribution<double> dist;

  Rand () : rg (rd ()), dist (0, 1) {};

  double next () {return dist (rg);}
};

struct Vec {
  double x;
  double y;

  Vec () : x (0), y (0) {};
  Vec (double _x, double _y) : x (_x), y (_y) {};
  Vec (const Vec& v) : x (v.x), y (v.y) {};
};

struct Particle {
  Vec p;
  Vec v;

  std::array<double, 1> data = {1.};

  Particle (Vec p) : p (p) {};
  Particle (Vec p, Vec v) : p (p), v (v) {};
  Particle (Vec p, Vec v, std::array<double, 1> data) : p (p), v (v), data (data) {};

  double operator[] (const size_t ID) const {return data [ID];};
  double& operator[] (const size_t ID) {return data [ID];};
};

struct id {
  int i;
  int j;

  id (int i, int j) : i (i), j (j) {};
};

Vec RK3 (const Vec& pos, const double dt, const VectorGrid& v);

#endif
