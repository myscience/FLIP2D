#include "include/utils.h"
#include "include/vector_grid.h"

Vec RK3 (const Vec& pos, const double dt, const VectorGrid& v) {
  double u1 = v.x.lerp (pos) / v.dx;
  double v1 = v.y.lerp (pos) / v.dy;

  double x1 = pos.x - 0.5 * dt * u1;
  double y1 = pos.y - 0.5 * dt * v1;

  double u2 = v.x.lerp (Vec (x1, y1)) / v.dx;
  double v2 = v.y.lerp (Vec (x1, y1)) / v.dy;

  double x2 = pos.x - 0.75 * dt * u2;
  double y2 = pos.y - 0.75 * dt * v2;

  double u3 = v.x.lerp (Vec (x2, y2)) / v.dx;
  double v3 = v.y.lerp (Vec (x2, y2)) / v.dy;

  Vec f_pos (pos.x, pos.y);

  f_pos.x += dt * (2. * u1 + 3. * u2 + 4. * u3) / 9.;
  f_pos.y += dt * (2. * v1 + 3. * v2 + 4. * v3) / 9.;

  return f_pos;
}
