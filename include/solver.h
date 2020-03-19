#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include <assert.h>
#include <math.h>
#include <tuple>

#include "grid.h"
#include "vector_grid.h"
#include "levelset.h"

struct Stencil {
  int w, h;

  double* diag = nullptr;
  double* plusX = nullptr;
  double* plusY = nullptr;

  const LevelSet& fls;

  Stencil (int w, int h, const LevelSet& ls) : w (w), h (h), diag (new double [w * h]()),
    plusX  (new double [w * h]()), plusY  (new double [w * h]()), fls (ls) {};

  ~Stencil () {
    if (diag   != nullptr) {delete[] diag;   diag   = nullptr;}
    if (plusX  != nullptr) {delete[] plusX;  plusX  = nullptr;}
    if (plusY  != nullptr) {delete[] plusY;  plusY  = nullptr;}
  }

  Grid operator* (Grid& u) const {
    Grid out (u.w, u.h);

    for (int j = 0; j < u.h; j++) {
      for (int i = 0; i < u.w; i++) {
        int _idx = i + j * w;

        if (fls [_idx] > 0) continue;

        id idx (i, j);

        double t = diag [_idx] * u[idx];

        if (fls [id (i + 1, j)] < 0) t += plusX [_idx] * u [id (i + 1, j)];
        if (fls [id (i, j + 1)] < 0) t += plusY [_idx] * u [id (i, j + 1)];
        if (fls [id (i - 1, j)] < 0) t += plusX [_idx - 1] * u [id (i - 1, j)];
        if (fls [id (i, j - 1)] < 0) t += plusY [_idx - w] * u [id (i, j - 1)];

        out [_idx] = t;
      }
    }

    return out;
  };
};

class Solver {
public:
  Solver (int w, int h, double dx, double dy, double dt, const LevelSet& ls);
  ~Solver ();

  VectorGrid project (const LevelSet& ls, const LevelSet& sls, VectorGrid& v);

public:
  Grid solveCG (const LevelSet& fls, const LevelSet& sls, const VectorGrid& v);
  Grid solvePCG (const LevelSet& fls, const LevelSet& sls, const VectorGrid& v);
  Grid CG (const Stencil& A, const Grid& rhs, double eps);
  Grid PCG (const Stencil& A, const Grid& rhs, const Grid& iE, const LevelSet& ls, double eps);
  Grid MIC (const Stencil& A, const LevelSet& ls, const double tau, const double sigma);
  Grid Psolve (const Grid& P, const Grid& r, const Stencil& A, const LevelSet& ls);

  Grid buildRhs (const Grid& u, double Vo, double mu);
  Grid buildCNrhs (const Grid& u, double delta);

  Grid buildRhoRhs (const Grid& rho, const VectorGrid& p);
  VectorGrid buildPRhs (const Grid& rho, const VectorGrid& p);

  /* ADIM Method */
  Grid buildXRhs (const Grid& u);
  Grid buildYRhs (const Grid& u);

  void buildStencil (const LevelSet& ls, const LevelSet& sls);

  void init ();
  double Vol (const Grid& u);

  int w;
  int h;
  int size;

  int maxIt = 1000;

  double dx;
  double dy;
  double dt;

  /* Solver parameters */
  double rho = 1.;

  Stencil A;
};

#endif
