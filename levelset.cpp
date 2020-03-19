#include "include/levelset.h"

LevelSet::LevelSet (const int w, const int h) : Grid (w, h) {};
LevelSet::LevelSet (std::string filename) : Grid (filename) {};
LevelSet::~LevelSet () {};

LevelSet& LevelSet::operator&= (const LevelSet& ls) {
  assert (w == ls.w and h == ls.h);

  for (int i = 0; i < size; i++) {
    if (ls [i] < 0) data [i] = 1.;
  }

  return *this;
}

void LevelSet::fit (const Particles& P) {
  /* Set current distance map phi to infinity and label to UNKWN */
  for (int i = 0; i < size; i++) {
    data  [i] = 1e20;
    _data [i] = UNKWN;
  }

  /* Here we initialize the array near the input geometry */
  size_t p_count = 0;
  for (const auto& p : P.Ps) {
    /* Here we locate particle p on the grid */
    int ix = (int) p.p.x;
    int iy = (int) p.p.y;

    double d = sqrt ((p.p.x - offX - ix) * (p.p.x - offX - ix) +
                     (p.p.y - offY - iy) * (p.p.y - offY - iy));

    /* Here we update current distance estimate */
    if (d < data [ix + iy * w]) {
      data  [ix + iy * w] = d;
      _data [ix + iy * w] = p_count;
    }

    p_count++;
  }

  /* Here we propagate closest point information to the rest of the grid */
  /* We make use of the fast sweep method: we iterate over our grid four
   * times: i ascending, j ascending, j descending, i descending. */
  /* I ascending */
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      /* We check idx neighbours */
      id _n (i, j + 1);
      id _s (i, j - 1);
      id _e (i + 1, j);
      id _w (i - 1, j);

      std::array<id, 4> neigh = {_n, _s, _w, _e};

      for (id idx : neigh) {
        size_t e = (*this) (idx);

        /* Check if we already know closest point of e */
        if (e >= (size_t) UNKWN) continue;

        Vec p = P.Ps.at (e).p;
        double d = sqrt ((p.x - offX - i) * (p.x - offX - i) +
                         (p.y - offY - j) * (p.y - offY - j));

        /* Check if update is needed */
        if (d < data [i + j * w]) {
          data  [i + j * w] = d;
          _data [i + j * w] = e;
        }
      }
    }
  }

  /* J ascending */
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      /* We check idx neighbours */
      id _n (i, j + 1);
      id _s (i, j - 1);
      id _e (i + 1, j);
      id _w (i - 1, j);

      std::array<id, 4> neigh = {_n, _s, _w, _e};

      for (id idx : neigh) {
        size_t e = (*this) (idx);

        /* Check if we already know closest point of e */
        if (e >= (size_t) UNKWN) continue;

        Vec p = P.Ps.at (e).p;
        double d = sqrt ((p.x - offX - i) * (p.x - offX - i) +
                         (p.y - offY - j) * (p.y - offY - j));

        /* Check if update is needed */
        if (d < data [i + j * w]) {
          data  [i + j * w] = d;
          _data [i + j * w] = e;
        }
      }
    }
  }

  /* J descending */
  for (int i = w - 1; i >= 0; i--) {
    for (int j = h - 1; j >= 0; j--) {
      /* We check idx neighbours */
      id _n (i, j + 1);
      id _s (i, j - 1);
      id _e (i + 1, j);
      id _w (i - 1, j);

      std::array<id, 4> neigh = {_n, _s, _w, _e};

      for (id idx : neigh) {
        size_t e = (*this) (idx);

        /* Check if we already know closest point of e */
        if (e >= (size_t) UNKWN) continue;

        Vec p = P.Ps.at (e).p;
        double d = sqrt ((p.x - offX - i) * (p.x - offX - i) +
                         (p.y - offY - j) * (p.y - offY - j));

        /* Check if update is needed */
        if (d < data [i + j * w]) {
          data  [i + j * w] = d;
          _data [i + j * w] = e;
        }
      }
    }
  }

  /* I descending */
  for (int j = h - 1; j >= 0; j--) {
    for (int i = w - 1; i >= 0; i--) {
      /* We check idx neighbours */
      id _n (i, j + 1);
      id _s (i, j - 1);
      id _e (i + 1, j);
      id _w (i - 1, j);

      std::array<id, 4> neigh = {_n, _s, _w, _e};

      for (id idx : neigh) {
        size_t e = (*this) (idx);

        /* Check if we already know closest point of e */
        if (e >= (size_t) UNKWN) continue;

        Vec p = P.Ps.at (e).p;
        double d = sqrt ((p.x - offX - i) * (p.x - offX - i) +
                         (p.y - offY - j) * (p.y - offY - j));

        /* Check if update is needed */
        if (d < data [i + j * w]) {
          data  [i + j * w] = d;
          _data [i + j * w] = e;
        }
      }
    }
  }

  /* Here we subtract the particle radius */
  (*this) -= P.r;

  /* Here we fill holes */
  fill_holes ();
}

void LevelSet::redistance () {

}

void LevelSet::fill_holes () {
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      id _n (i, j + 1);
      id _s (i, j - 1);
      id _e (i + 1, j);
      id _w (i - 1, j);

      std::array<id, 4> neigh = {_n, _s, _w, _e};

      /* Compute the average of the neighbours */
      double old = data [i + j * w];
      double avg = 0.;
      size_t N = 0;

      for (const auto& idx : neigh)
        if ((*this) (idx) < UNKWN) {avg += (*this) [idx]; N++;}

      assert (N != 0);
      avg /= N;

      /* Check if hole filling is needed */
      data [i + j * w] = avg < old ? avg : old;
    }
  }
}

/* FIXME: This routine is probably brocken */
Vec LevelSet::find_closest (const Vec &x) const {
  /* Only move if inside */
  if (lerp (x) > 0) return x;

  const int N = 20;
  const int M = 50;

  Vec p (x);
  Vec d = lerpGrad (p);
  double phi_p = lerp (p);

  for (int i = 0; i < N; i++) {
    double alpha = 1.;

    for (int j = 0; j < M; j++) {
      Vec q (p.x - alpha * phi_p * d.x, p.y - alpha * phi_p * d.y);

      if (abs (lerp (q)) < abs (phi_p)) {
        p.x = q.x; p.y = q.y;

        phi_p = lerp (q);
        d = lerpGrad (q);

        if (abs (phi_p) < 1e-3) return p;
      } else {
        alpha *= 0.7;
      }
    }
  }

  return Vec (-1, -1);
}

Grid LevelSet::volume () const {
  Grid V (w, h);

  for (int j = 0; j < h - 1; j++) {
    for (int i = 0; i < w - 1; i++) {
      int idx = i + j * w;

      double phi1 = data [idx], phi2 = data [idx + 1],
             phi3 = data [idx + w], phi4 = data [idx + 1 + w];

      double phi0 = 0.5 * (phi1 + phi2 + phi3 + phi4);

      double v = _volume (phi0, phi1, phi2, phi3, phi4);
      V [idx] = v < 0.01 ? 0 : v;
    }
  }

  return V;
}

double LevelSet::_volume (double phi0, double phi1, double phi2, double phi3, double phi4) const {
  return 0.25 * (_volume (phi0, phi1, phi2) + _volume (phi0, phi2, phi3) +
                 _volume (phi0, phi3, phi4) + _volume (phi0, phi4, phi1));
}

double LevelSet::_volume (double a, double b, double c) const {
  std::array<double, 3> V = {a, b, c};
  std::sort (std::begin (V), std::end (V));

  if (V [0] >= 0) return 0;
  if (V [2] <= 0) return 1;
  if (V [1] <= 0 and V [2] > 0) return (V [2] / (V [2] - V [0])) * (V [2] / (V [2] - V [1]));
  if (V [0] <= 0 and V [1] > 0) return 1. - (V [0] / (V [0] - V [1])) * (V [0] / (V [0] - V [2]));

  return -1;
}
