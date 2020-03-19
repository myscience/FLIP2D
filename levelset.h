#ifndef LEVELSET_H
#define LEVELSET_H

#include <math.h>
#include <array>
#include <algorithm>

#include "grid.h"
#include "particles.h"

/* This is our LevelSet class that extends the Grid class */
class LevelSet : public Grid {
public:
  LevelSet (const int w, const int h);
  LevelSet (std::string filename);
  ~LevelSet ();

  void fit (const Particles& P);
  void redistance ();
  void fill_holes ();

  Grid volume () const;

  LevelSet& operator&= (const LevelSet& ls);

  Vec find_closest (const Vec& x) const;

private:
  double _volume (double a, double b, double c, double d, double e) const;
  double _volume (double a, double b, double c) const;
};

#endif
