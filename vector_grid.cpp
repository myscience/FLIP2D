#include "include/vector_grid.h"

VectorGrid::VectorGrid (int w, int h, double dx, double dy) :
  w (w), h (h), dx (dx), dy (dy), x (w + 1, h, 0., 0.5), y (w, h + 1, 0.5, 0.) {};

VectorGrid::VectorGrid (double dx, double dy, Grid& _x, Grid& _y) :
  w (_x.w), h (_x.h), dx (dx), dy (dy), x (_x), y (_y) {
    assert (_x.w == _y.w && _x.h == _y.h);
  };

VectorGrid::VectorGrid (int w, int h, double dx, double dy, std::string init_x, std::string init_y) :
  w (w), h (h), dx (dx), dy (dy), x (0., 0.5, init_x), y (0.5, 0., init_y) {
  assert (x.w == y.w + 1 && x.h + 1 == y.h);
};

VectorGrid::VectorGrid (const VectorGrid& v) : w (v.w), h (v.h),
  dx (v.dx), dy (v.dy), x (v.x), y (v.y) {};


VectorGrid::~VectorGrid () {};

Grid VectorGrid::operator*(VectorGrid &b) const {
  Grid out (w, h);

  for (int i = 0; i < w * h; i++)
    out [i] = x [i] * b.x [i] + y [i] * b.y [i];

  return out;
}

Grid VectorGrid::operator*(const VectorGrid &b) const {
  Grid out (w, h);

  for (int i = 0; i < w * h; i++)
    out [i] = x [i] * b.x [i] + y [i] * b.y [i];

  return out;
}

VectorGrid VectorGrid::operator* (const double c) const {
  VectorGrid out (*this);

  out.x *= c;
  out.y *= c;

  return out;
}

VectorGrid VectorGrid::operator- (const VectorGrid& v) const {
  VectorGrid out (*this);

  out.x -= v.x;
  out.y -= v.y;

  return out;
}

VectorGrid& VectorGrid::operator+= (const VectorGrid v) {
  assert (dx == v.dx && dy == v.dy);

  x += v.x;
  y += v.y;

  return *this;
}

VectorGrid& VectorGrid::operator= (const VectorGrid& v) {
  assert (dx == v.dx && dy == v.dy);

  x = v.x;
  y = v.y;

  return *this;
}

VectorGrid VectorGrid::grad (const Grid& a) const {
  VectorGrid grad (a.w, a.h, dx, dy);

  double idx = 1. / dx;
  double idy = 1. / dy;

  for (int i = 0; i < a.w; i++) {
    for (int j = 0; j < a.h; j++) {
      grad.x [id (i, j)] = idx * (a [id (i + 1, j)] - a [id (i, j)]);
      grad.y [id (i, j)] = idy * (a [id (i, j + 1)] - a [id (i, j)]);
    }
  }

  return grad;
}

Grid VectorGrid::ndiv (const LevelSet& fls, const LevelSet& sls) const {
  Grid ndiv (w, h);

  double ivdx = 1. / dx;
  double ivdy = 1. / dy;

  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      int idx = i + j * w;

      if (fls [idx] > 0) continue;

      ndiv [idx] = -ivdx * (x [id (i + 1, j)] - x [id (i, j)])
                   -ivdy * (y [id (i, j + 1)] - y [id (i, j)]);

      /* Modify RHS to account for solid velocities */
      /* NOTE: We assume that solid has no velocities */
      if (sls [idx - 1] < 0) ndiv [idx] -= ivdx * (x [id (i, j)] - 0.);
      if (sls [idx + 1] < 0) ndiv [idx] += ivdx * (x [id (i + 1, j)] - 0.);
      if (sls [idx - w] < 0) ndiv [idx] -= ivdy * (y [id (i, j)] - 0.);
      if (sls [idx + w] < 0) ndiv [idx] += ivdx * (y [id (i, j + 1)] - 0.);
    }
  }

  return ndiv;
}

double VectorGrid::norm_2 () const {
  return x.norm_2 () + y.norm_2 ();
}

void VectorGrid::advect (const double dt, VectorGrid &v) {
  x.advect (dt, v);
  y.advect (dt, v);
}

void VectorGrid::swap(VectorGrid& new_v) {
  using std::swap;

  swap(*this, new_v);
}

void VectorGrid::fit (const Particles& P) {
  /* First we clear the grid */
  x.clear ();
  y.clear ();

  /* This is the weight grid */
  Grid Wx (w + 1, h);
  Grid Wy (w, h + 1);

  /* Here we loop over the particles */
  for (const auto& P : P.Ps) {
    /* We extract particles position */
    Vec p = P.p;
    Vec v = P.v;

    id idx ((int) (p.x), (int) (p.y));
    id _n (idx.i, idx.j + 1);
    id _s (idx.i, idx.j - 1);
    id _e (idx.i + 1, idx.j);
    id _w (idx.i - 1, idx.j);

    id _nn (idx.i, idx.j + 2);
    id _ss (idx.i, idx.j - 2);
    id _ee (idx.i + 2, idx.j);
    id _ww (idx.i - 2, idx.j);

    id _ne (idx.i + 1, idx.j + 1);
    id _nw (idx.i - 1, idx.j + 1);
    id _se (idx.i + 1, idx.j - 1);
    id _sw (idx.i - 1, idx.j - 1);

    // std::array<id, 5> neigh = {idx, _n, _s, _w, _e};
    std::array<id, 13> neigh = {idx, _n, _s, _w, _e, _nn, _ss, _ww, _ee, _ne, _nw, _se, _sw};

    for (const auto& idx : neigh) {
      /* Avoid out-of-border idxs */
      if (x (idx) == x.nulllabel) continue;

      /* Add the contribution */
      double k = x.Bspline (P.p, idx);
      x [idx] += v.x * k;

      Wx [idx] += k;
    }

    for (const auto& idx : neigh) {
      /* Avoid out-of-border idxs */
      if (y (idx) == y.nulllabel) continue;

      /* Add the contribution */
      double k = y.Bspline (P.p, idx);
      y [idx] += v.y * k;

      Wy [idx] += k;
    }
  }

  /* Here we divide by the weights */
  for (int i = 0; i < x.size; i++) {
    if (Wx [i] != 0.) {
      x [i] /= Wx [i];
      x (i) = KNOWN;
    } else {
      x (i) = UNKWN;
    }
  }

  for (int i = 0; i < y.size; i++) {
    if (Wy [i] != 0.) {
      y [i] /= Wy [i];
      y (i) = KNOWN;
    } else {
      y (i) = UNKWN;
    }
  }

  /* Finally, we extrapolate the inferred values */
  x.extrapolate ();
  y.extrapolate ();
}
