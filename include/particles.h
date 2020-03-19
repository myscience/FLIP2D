#ifndef PARTICLES_H
#define PARTICLES_H

#include "utils.h"

class LevelSet;

class Particles {
public:
  Particles (const double r, Grid& count, std::string filename);
  ~Particles ();

  void advect (const VectorGrid& v, const LevelSet& sls, const double dt);
  void fit (const VectorGrid& v_old, const VectorGrid& v_new, double alpha);
  void fit (const Grid& q_old, const Grid& q_new, size_t ID, double alpha);
  void adjust (const VectorGrid& v, const Grid& rho, const LevelSet& fls, const LevelSet& sls);

  void recount () const;
  void prune (const LevelSet& fls);
  void seed (const LevelSet& fls, const LevelSet& sls, const VectorGrid& v, const Grid& rho);

  int w;
  int h;
  double r;

  Grid& count;
  Rand R;

  std::vector<Particle> Ps;

  const size_t minPerCell = 3;
  const size_t avgPerCell = 8;
  const size_t maxPerCell = 12;
};

#endif
