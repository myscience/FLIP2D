#include "include/grid.h"
#include "include/levelset.h"
#include "include/particles.h"

Grid::Grid (const int& w, const int& h) : w (w), h (h), size (w * h),
  offX (0.5), offY (0.5), data (new double[size]()), _data (new uint[size]()) {};

Grid::Grid (const int& w, const int& h, const int offX, const int offY) : w (w), h (h), size (w * h),
  offX (offX), offY (offY), data (new double[size]()), _data (new uint[size]()) {};

Grid::Grid (std::string filename) : offX (0.5), offY (0.5) {
  std::ifstream in (filename);

  in >> *this;
};

Grid::Grid (const double dx, const double dy, std::string filename) : offX (dx), offY (dy) {
  std::ifstream in (filename);

  in >> *this;
};

Grid::Grid (const Grid& a) : w (a.w), h (a.h), size (w * h), offX (a.offX), offY (a.offY),
data (new double[size]()), _data (new uint[size]()) {
  std::copy (a.data, a.data + a.size, data);
};

Grid::~Grid () {
  if (data != nullptr) {
    delete[] data;
    data = nullptr;
  }

  if (_data != nullptr) {
    delete[] _data;
    _data = nullptr;
  }
};

double& Grid::operator[] (int i) {
  /* We check boundary condition */
  if (i >= 0 && i < size) return data [i];

  return nullvalue;
}

const double& Grid::operator[] (int i) const {
  /* We check boundary condition */
  if (i >= 0 && i < size) return data [i];

  return nullvalue;
}

uint& Grid::operator() (int i) {
  /* We check boundary condition */
  if (i >= 0 && i < size) return _data [i];

  return nulllabel;
}

const uint& Grid::operator() (int i) const {
  /* We check boundary condition */
  if (i >= 0 && i < size) return _data [i];

  return nulllabel;
}

double& Grid::operator[] (id idx) {
  /* We check boundary condition */
  if ((idx.i > -1 && idx.i < w) &&
      (idx.j > -1 && idx.j < h)) return data [idx.j * w + idx.i];

  return nullvalue;
};

const double& Grid::operator[] (id idx) const {
  /* We check boundary condition */
  if ((idx.i > -1 && idx.i < w) &&
      (idx.j > -1 && idx.j < h)) return data [idx.j * w + idx.i];

  return nullvalue;
};


uint& Grid::operator() (id idx) {
  /* We check boundary condition */
  if ((idx.i > -1 && idx.i < w) &&
      (idx.j > -1 && idx.j < h)) return _data [idx.j * w + idx.i];

  return nulllabel;
}

const uint& Grid::operator() (id idx) const {
  /* We check boundary condition */
  if ((idx.i > -1 && idx.i < w) &&
      (idx.j > -1 && idx.j < h)) return _data [idx.j * w + idx.i];

  return nulllabel;
}

Grid& Grid::operator= (const Grid& b) {
  assert (size = b.size);

  std::copy (b.data, b.data + b.size, data);
  std::copy (b._data, b._data + b.size, _data);

  return *this;
}

/* This is the Grid equivalent of the dot-product */
double Grid::operator* (const Grid& b) const {
  assert (size == b.size);
  double r = 0;

  for (int i = 0; i < size; i++)
    r += data[i] * b [i];

  return r;
}

/* This is multiplication with scalar */
Grid Grid::operator* (const double c) {
  Grid out (*this);

  for (int i = 0; i < size; i++)
    out [i] *= c;

  return out;
}

Grid Grid::operator* (const double c) const {
  Grid out (*this);

  for (int i = 0; i < size; i++)
    out [i] *= c;

  return out;
}

Grid Grid::operator- (const double b) const {
  Grid out (*this);

  for (int i = 0; i < size; i++)
    out [i] -= b;

  return out;
}


/* This is multiplication with scalar */
Grid& Grid::operator*= (const double c) {
  for (int i = 0; i < size; i++)
    data [i] *= c;

  return *this;
}

Grid Grid::operator+ (const Grid& b) const {
  assert (size == b.size);

  Grid out (b);
  for (int i = 0; i < size; i++)
    out [i] += data [i];

  return out;
}

Grid Grid::operator- (const Grid& b) const {
  assert (size == b.size);

  Grid out (b);

  for (int i = 0; i < size; i++)
    out [i] -= data [i];

  return out;
}

Grid& Grid::operator+= (const double c) {
  for (int i = 0; i < size; i++)
    data [i] += c;

  return *this;
}

Grid& Grid::operator+= (const Grid& b) {
  assert (size == b.size);

  for (int i = 0; i < size; i++)
    data [i] += b.data [i];

  return *this;
}

Grid& Grid::operator-= (const double c) {
  for (int i = 0; i < size; i++)
    data [i] -= c;

  return *this;
}

Grid& Grid::operator-= (const Grid& b) {
  assert (size == b.size);

  for (int i = 0; i < size; i++)
    data [i] -= b.data [i];

  return *this;
}

bool Grid::operator== (double c) const {
  bool ans = true;

  for (int i = 0; i < w * h; i++)
    ans &= data [i] == c;

  return ans;
}

double Grid::norm_inf () const {
  double max = -1e50;

  for (int i = 0; i < size; i++)
    if (fabs (data [i]) > max) max = fabs (data [i]);

  return max;
}

double Grid::norm_2 () const {
  double norm = 0;

  for (int i = 0; i < size; i++)
    norm += data [i] * data [i];

  return norm / size;
}

Grid Grid::smooth (int sx, int sy) {
  assert (w % sx <= 1 && h % sy <= 1);

  int nw = w / sx;
  int nh = h / sy;

  Grid out (nw, nh);

  for (int k = 0; k < nw; k++) {
    for (int l = 0; l < nh; l++) {
      double mean = 0.;

      for (int i = k * sx; i < (k + 1) * sx; i++)
        for (int j = l * sy; j < (l + 1) * sy; j++)
          mean += data [i + j * w];

      out [id (k, l)] = mean / (sx * sy);
    }
  }

  return out;
}

double Grid::lerp (const Vec& pos) const {
  using std::max;
  using std::min;

  double x = min (max (pos.x - offX, 0.), w - 1.0001);
  double y = min (max (pos.y - offY, 0.), h - 1.0001);

  int ix = (int) x;
  int iy = (int) y;
  x -= ix;
  y -= iy;

  /* Here we collect the neighbouring value, checking for borders */
  int ix_ = ix + 1;
  int iy_ = iy + 1;

  double x00 = data [ix + iy  * w], x10 = data [ix_ + iy  * w];
  double x01 = data [ix + iy_ * w], x11 = data [ix_ + iy_ * w];

  return _lerp (_lerp (x00, x10, x), _lerp (x01, x11, x), y);
}

Vec Grid::lerpGrad (const Vec& pos) const {
  using std::max;
  using std::min;

  double x = min (max (pos.x - offX, 0.), w - 1.0001);
  double y = min (max (pos.y - offY, 0.), h - 1.0001);

  int ix = (int) x;
  int iy = (int) y;
  x -= ix;
  y -= iy;

  /* Here we collect the neighbouring value, checking for borders */
  int ix_ = ix + 1;
  int iy_ = iy + 1;

  double x00 = data [ix + iy  * w], x10 = data [ix_ + iy  * w];
  double x01 = data [ix + iy_ * w], x11 = data [ix_ + iy_ * w];

  /* Here we compute the Gradient */
  double gx0 = (x10 - x00), gx1 = (x11 - x01);
  double gy0 = (x01 - x00), gy1 = (x11 - x10);

  return Vec (_lerp (gx0, gx1, y), _lerp (gy0, gy1, x));
}

double Grid::cerp (const Vec& pos) const {
  using std::max;
  using std::min;

  int ix = (int) pos.x - offX;
  int iy = (int) pos.y - offY;

  double x = pos.x - offX - ix;
  double y = pos.y - offY - iy;

  /* Here we collect the id for the four neighbours */
  int x0 = max (ix - 1, 0), x1 = ix, x2 = min (ix + 1, w - 1), x3 = min (ix + 2, w - 1);
  int y0 = max (iy - 1, 0), y1 = iy, y2 = min (iy + 1, h - 1), y3 = min (iy + 2, h - 1);

  double q0 = _cerp (data [x0 + y0 * w], data [x1 + y0 * w], data [x2 + y0 * w], data [x3 + y0 * w], x);
  double q1 = _cerp (data [x0 + y1 * w], data [x1 + y1 * w], data [x2 + y1 * w], data [x3 + y1 * w], x);
  double q2 = _cerp (data [x0 + y2 * w], data [x1 + y2 * w], data [x2 + y2 * w], data [x3 + y2 * w], x);
  double q3 = _cerp (data [x0 + y3 * w], data [x1 + y3 * w], data [x2 + y3 * w], data [x3 + y3 * w], x);

  return _cerp (q0, q1, q2, q3, y);
}

double Grid::_cerp (double a, double b, double c, double d, double x) const {
  using std::max;
  using std::min;

  double xsq = x * x;
  double xcu = xsq * x;

  double minV = min (a, min (b, min (c, d)));
  double maxV = max (a, max (b, max (c, d)));

  double t = a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu) +
             b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu) +
             c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu) +
             d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

  return min (max (t, minV), maxV);
}

double Grid::_lerp (double a, double b, double x) const {
  return a * (1. - x) + b * x;
}

double Grid::_h (double r) const {
  if (r >= - 1.5 && r < -0.5) return 0.5 * (r + 1.5) * (r + 1.5);
  else if (r >= -0.5 && r < 0.5) return 0.75 - (r * r);
  else if (r >= 0.5 && r < 1.5) return 0.5 * (1.5 - r) * (1.5 - r);
  else return 0.;
}

double Grid::Bspline (const Vec& v, const id& idx) {
  /* Here we compute rx and ry */
  double rx = (v.x - offX) - idx.i;
  double ry = (v.y - offY) - idx.j;

  return _h (rx) * _h (ry);
}

void Grid::clear () {
  for (int i = 0; i < size; i++) {
    data  [i] = 0;
    _data [i] = UNKWN;
  }
}

void Grid::advect (const double dt, const VectorGrid& v) {
  using std::swap;

  Grid new_grid (w, h);

  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      Vec pos = RK3 (Vec (i, j), -dt, v);

      new_grid [id (i, j)] = cerp (pos);
    }
  }

  // *this = new_grid;
  swap (*this, new_grid);
}

void Grid::fit (const Particles& P, size_t ID) {
  /* First we clear the grid */
  clear ();

  /* This is the weight grid */
  Grid W (w, h);

  /* Here we loop over the particles */
  for (const auto& P : P.Ps) {
    /* We extract particles position */
    Vec p = P.p;

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

    for (const auto& _idx : neigh) {
      /* Avoid out-of-border idxs */
      if ((*this) (_idx) == nulllabel) continue;

      /* Add the contribution */
      double k = Bspline (P.p, _idx);
      (*this) [_idx] += P [ID] * k;
      // W [_idx] += 1;
    }
  }

  /* Here we divide by the weights */
  // P.recount ();

  for (int i = 0; i < size; i++) {
    int count = P.count [i];

    if (count != 0.) {
      data [i] /= 1;
      _data [i] = KNOWN;
    } else {
      _data [i] = UNKWN;
    }
  }

  // extrapola  te ();
}

/* FIXME: THIS FUNCTION DOES NOT WORK! */
void Grid::extrapolate () {
  /* Here we allocate the Wavefront vector W */
  std::vector<id> W;

  /* First we init the wavefront */
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      id _n (i, j + 1);
      id _s (i, j - 1);
      id _e (i + 1, j);
      id _w (i - 1, j);

      std::array<id, 4> neigh = {_n, _s, _w, _e};

      bool flag = false;
      for (const auto& idx : neigh) flag |= (*this) (idx) == KNOWN;

      /* Check if (i, j) belongs to wavefront */
      if (_data [i + j * w] == UNKWN && flag) {
        _data [i + j * w] = 1;
        W.push_back (id (i, j));
      }
    }
  }

  size_t t = 0;
  while (t < W.size ()) {
    id idx = W [t];

    id _n (idx.i, idx.j + 1);
    id _s (idx.i, idx.j - 1);
    id _e (idx.i + 1, idx.j);
    id _w (idx.i - 1, idx.j);

    std::array<id, 4> neigh = {_n, _s, _w, _e};

    size_t N = 0;
    double F = 0;

    for (const auto& _idx : neigh) {
      if ((*this) (_idx) < (*this) (idx)) {
        F += (*this) [_idx];
        N++;
      }

      /* If _idx was previously UNKWN, now it can be added to the wavefront */
      if ((*this) (_idx) == UNKWN) {
        (*this) (_idx) = (*this) (idx) + 1;
        W.push_back (_idx);
      }
    }

    /* Set this grid idx to the average of reliable neighbours */
    assert (N != 0);
    (*this) [idx] = F / N;

    t++;
  }
}

std::ostream& operator <<(std::ostream& stream, const Grid& grid) {
  for (int j = 0; j < grid.h; j++)
    for (int i = 0; i < grid.w; i++)
      stream << grid.data [i + j * grid.w] << ' ';// << ((i == grid.w - 1) ? '\n' : ' ');

  stream << '\n';

  return stream;
};

std::istream& operator >>(std::ifstream& stream, Grid& grid) {
  /* Here we parse the file */
  char c;
  stream >> c >> grid.w >> grid.h;

  grid.size = grid.w * grid.h;

  grid.data = new double[grid.size]();
  grid._data = new uint [grid.size]();

  for (int j = 0; j < grid.h; j++) {
    for (int i = 0; i < grid.w; i++) {
      int idx = i + j * grid.w;

      stream >> grid.data [idx];
      grid._data [idx] = grid.data [idx] > 0 ? FLUID : EMPTY;
    }
  }

  return stream;
};

unsigned char* operator<< (unsigned char* rgba, const Grid& grid) {
  for (int j = 0; j < grid.h; j++) {
    for (int i = 0; i < grid.w; i++) {
      int idx_rgba = 4 * (i + j * grid.w);

      double shade = 1. - grid [i + j * grid.w];

      shade = std::min (std::max (shade, 0.0), 1.0);
      rgba [idx_rgba + 0] = (int) (shade * 255.0);
      rgba [idx_rgba + 1] = (int) (shade * 255.0);
      rgba [idx_rgba + 2] = (int) (shade * 255.0);
      rgba [idx_rgba + 3] = 0xFF;
    }
  }

  return rgba;
}
