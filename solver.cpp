#include "include/solver.h"

Solver::Solver (int w, int h, double dx, double dy, double dt, const LevelSet& ls) :
  w (w), h (h), size (w * h), dx (dx), dy (dy), dt (dt), A (w, h, ls) {};

Solver::~Solver() {};

void Solver::buildStencil (const LevelSet& ls, const LevelSet& sls) {
  double scaleX = dt / (rho * 2 * dx * dx);
  double scaleY = dt / (rho * 2 * dy * dy);

  double scale = scaleX + scaleY;

  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      int idx = i + j * w;

      /* Here we clear previously computed values */
      A.diag [idx] = A.plusX [idx] = A.plusY [idx] = 0.;

      if (ls [idx] > 0) continue;

      id l (i - 1, j);
      id r (i + 1, j);
      id b (i, j - 1);
      id t (i, j + 1);

      double ilidx = ls [idx] > -1e-6 ? -1e6 : 1. / ls [idx];

      /* This is the x direction */
      if (ls [l] < 0)  A.diag [idx] += scale;
      else if (ls [l] > 0 && sls [l] > 0) A.diag [idx] += scale * (ls [idx] - ls [l]) * ilidx;

      if (ls [r] < 0) {A.diag [idx] += scale; A.plusX [idx] = -scale;}
      else if (ls [r] > 0 && sls [r] > 0) A.diag [idx] += scale * (ls [idx] - ls [r]) * ilidx;

      /* This is the y direction */
      if (ls [b] < 0)  A.diag [idx] += scale;
      else if (ls [b] > 0 && sls [b] > 0) A.diag [idx] += scale * (ls [idx] - ls [b]) * ilidx;

      if (ls [t] < 0) {A.diag [idx] += scale; A.plusY [idx] = -scale;}
      else if (ls [t] > 0 && sls [t] > 0) A.diag [idx] += scale * (ls [idx] - ls [t]) * ilidx;
    }
  }
}

double Solver::Vol (const Grid& u) {
  double V = 0;
  for (int i = 0; i < size; i++)
      V += u [i];

  return V / size;
}

Grid Solver::MIC (const Stencil& A, const LevelSet& ls, const double tau, const double sigma) {
  /* Here we build the Modified Incomplete Cholesky Preconditioner */
  /* We actually compute only the (inverse) of the diagonal of L, because for
   * our matrix A it can be shown that L = F E^-1 + E, where F is the low
   * triangular part of A. Thus we only compute E^-1 (for efficiency). */

   Grid iE (w, h);

   for (int j = 0; j < h; j++) {
     for (int i = 0; i < w; i++) {
       int idx = i + j * w;

       if (ls [idx] > 0) continue;

       double e = A.diag [idx];

       if (ls [id (i - 1, j)] < 0) {
         double px = A.plusX [idx - 1] * iE [idx - 1];
         double py = A.plusY [idx - 1] * iE [idx - 1];

         e -= (px * px + tau * px * py);
       }

       if (ls [id (i, j - 1)] < 0) {
         double px = A.plusX [idx - w] * iE [idx - w];
         double py = A.plusY [idx - w] * iE [idx - w];

         e -= (py * py + tau * px * py);
       }

       if (e < sigma * A.diag [idx]) e = A.diag [idx];

       iE [id (i, j)] = e > 1e-9 ? 1. / sqrt(e) : 0;
     }
   }

   return iE;
}

Grid Solver::Psolve (const Grid& P, const Grid& r, const Stencil& A, const LevelSet& ls) {
  /* First we solve Lq = r */
  Grid z (w, h);

  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      int idx = i + j * w;

      if (ls [idx] > 0) continue;


      double t = r [idx];

      if (ls [id (i - 1, j)] < 0)
        t -= A.plusX [idx - 1] * P [idx - 1] * z [idx - 1];

      if (ls [id (i, j - 1)] < 0)
        t -= A.plusY [idx - w] * P [idx - w] * z [idx - w];

      z [idx] = t * P [idx];
    }
  }

  /* Next solve L^T z = q */
  for (int j = h - 1; j >= 0; j--) {
    for (int i = w - 1; i >= 0; i--) {
      int idx = i + j * w;

      if (ls [idx] > 0) continue;

      double t = z [idx];

      if (ls [id (i + 1, j)] < 0)
        t -= A.plusX [idx] * P [idx] * z [idx + 1];

      if (ls [id (i, j + 1)] < 0)
        t -= A.plusY [idx] * P [idx] * z [idx + w];

      z [idx] = t * P [idx];
    }
  }

  return z;
}

Grid Solver::CG (const Stencil& A, const Grid& rhs, double eps) {
  /* This is our solution vector */
  Grid x (w, h);

  /* We init the residual vector */
  Grid r (rhs);

  /* Check if r is already close to zero */
  if (r == 0) return x;

  /* We init our guess to the rhs */
  Grid p (r);

  double sigma = r * r;

  for (int k = 0; k < maxIt; k++) {
    Grid z = A * p;

    double alpha = sigma / (p * z);

    x += p * alpha;
    r -= z * alpha;

    if (r.norm_inf () < eps) return x;

    double sigma_new = r * r;
    double beta = sigma_new / sigma;

    p = r + p * beta;
    sigma = sigma_new;
  }

  std::cerr << "WARNING: Exiting CG with Err: " << r.norm_inf () << '\n';

  return x;
};

Grid Solver::PCG (const Stencil& A, const Grid& rhs, const Grid& iE, const LevelSet& fls, double eps) {
  /* This is our solution vector */
  Grid p (w, h);

  /* We init the residual vector */
  Grid r (rhs);

  if (r == 0) return p;

  /* We init our guess to rhs */
  Grid z = Psolve (iE, r, A, fls);
  Grid s (z);

  double sigma = z * r;

  for (int k = 0; k < maxIt; k++) {
    Grid z = A * s;

    double alpha = sigma / (z * s);

    p += s * alpha;
    r -= z * alpha;

    if (r.norm_inf () < eps) return p;

    z = Psolve (iE, r, A, fls);

    double sigma_new = z * r;
    double beta = sigma_new / sigma;

    s = z + s * beta;
    sigma = sigma_new;
  }

  std::cerr << "WARNING: Exiting PCG with Err: " << r.norm_inf () << '\n';

  return p;
}

Grid Solver::solvePCG (const LevelSet& ls, const LevelSet& sls, const VectorGrid& v) {
  /* Here we combine the implicit CN method with the explicit RK4 method */
  buildStencil (ls, sls);

  /* Then we compute the negative divergence which is our RHS */
  Grid ndiv = v.ndiv (ls, sls);

  /* Then we build our Modified Incomplete Cholesky Preconditioner */
  Grid P = MIC (A, ls, 0.97, 0.25);

  /* Then we solve using the Conjugate Gradient */
  Grid p = PCG (A, ndiv, P, ls, 1e-4);

  /* Check if nan */
  assert (p [12] == p [12]);

  return p;
}

Grid Solver::solveCG (const LevelSet& fls, const LevelSet& sls, const VectorGrid& v) {
  /* First we build the Stencil of the linear system */
  buildStencil (fls, sls);

  /* Then we compute the negative divergence which is our RHS */
  Grid ndiv = v.ndiv (fls, sls);

  /* Then we solve using the Conjugate Gradient */
  Grid p = CG (A, ndiv, 1e-4);

  return p;
}

VectorGrid Solver::project (const LevelSet& fls, const LevelSet& sls, VectorGrid& v) {
  VectorGrid v_new (v);

  Grid p = solvePCG (fls, sls, v);

  double scaleX = dt / (rho * dx);
  double scaleY = dt / (rho * dy);

  for (int j = 0; j < h + 1; j++) {
    for (int i = 0; i < w + 1; i++) {
      id idx  = id (i, j);
      id id_l = id (i - 1, j);
      id id_b = id (i, j - 1);

      /* Here we update the u component of velocity */
      if (fls [id_l] < 0 or fls [idx] < 0) {
        v_new.x [idx] -= scaleX * (p [idx] - p [id_l]);

        /* Here we enforce boundary conditions */
        if (sls [idx] < 0) v_new.x [idx] = v_new.x [idx] > 0. ? 0. : v_new.x [idx];
        if (sls [id_l] < 0) v_new.x [idx] = v_new.x [idx] < 0. ? 0. : v_new.x [idx];

        /* Mark this grid position as known */
        v_new.x (idx) = KNOWN;
      } else {
        v_new.x (idx) = UNKWN;
      }

      /* Here we update the v component of velocity */
      if (fls [idx] < 0 or fls [id_b] < 0) {
        v_new.y [idx] -= scaleY * (p [idx] - p [id_b]);

        /* Here we enforce boundary conditions */
        if (sls [idx] < 0) v_new.y [idx] = v_new.y [idx] > 0. ? 0. : v_new.y [idx];
        if (sls [id_b] < 0) v_new.y [idx] = v_new.y [idx] < 0. ? 0. : v_new.y [idx];

        /* Mark this grid position as known */
        v_new.y (idx) = KNOWN;
      } else {
        v_new.y (idx) = UNKWN;
      }
    }
  }

  v_new.x.extrapolate ();
  v_new.y.extrapolate ();

  return v_new;
}
