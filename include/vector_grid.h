#ifndef VECTOR_GRID_H
#define VECTOR_GRID_H

#include "grid.h"
#include "levelset.h"

class VectorGrid {
public:
  VectorGrid (int w, int h, double dx, double dy);
  VectorGrid (double dx, double dy, Grid& x, Grid &y);
  VectorGrid (int w, int h, double dx, double dy, std::string in_x, std::string in_y);
  VectorGrid (const VectorGrid& v);

  ~VectorGrid ();

  /* Assignment operator */
  VectorGrid& operator= (const VectorGrid& v);

  Grid operator*(VectorGrid& b) const;
  Grid operator*(const VectorGrid& b) const;

  VectorGrid& operator+= (const VectorGrid v);

  VectorGrid operator* (const double dt) const;
  VectorGrid operator- (const VectorGrid& v) const;

  VectorGrid grad (const Grid& a) const;

  Grid ndiv (const LevelSet& fls, const LevelSet& sls) const;

  void advect (const double dt, VectorGrid& v);
  void fit (const Particles& P);
  double norm_2 () const;

  void swap (VectorGrid& v);

  int w;
  int h;
  double dx;
  double dy;

  Grid x;
  Grid y;
private:
};

#endif
