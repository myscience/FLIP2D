#include "include/particles.h"
#include "include/grid.h"
#include "include/levelset.h"
#include "include/vector_grid.h"

Particles::Particles (const double r, Grid& count, std::string filename) :
  r (r), count (count), R () {
  /* Here we parse the file with the initial fluid position to populate our
   * Particle vector */
  std::ifstream in (filename);

  char c;

  /* We parse width and height information */
  in >> c >> w >> h;
  assert (w == count.w && h == count.h);

  /* Here we reserve memory for the particles */
  Ps.reserve (w * h * avgPerCell);

  double value = 0.;
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      in >> value;

      if (value > 0.) continue;

      for (size_t c = 0; c < avgPerCell; c++) {
        /* Create new position for particles */
        double x = i + R.next ();
        double y = j + R.next ();

        Ps.push_back (Particle (Vec (x, y), Vec (0., 0)));
      }

      /* Set count to avgPerCell */
      count [id (i, j)] = avgPerCell;
    }
  }
};

Particles::~Particles () {};

void Particles::fit (const VectorGrid& v_old, const VectorGrid& v_new, double alpha) {
  VectorGrid delta = v_new - v_old;

  for (auto& P : Ps) {
    /* This is the PIC/FLIP update */
    P.v.x = alpha * v_new.x.lerp (P.p) + (1. - alpha) * (P.v.x + delta.x.lerp (P.p));
    P.v.y = alpha * v_new.y.lerp (P.p) + (1. - alpha) * (P.v.y + delta.y.lerp (P.p));
  }
}


void Particles::fit (const Grid& q_old, const Grid& q_new, size_t ID, double alpha) {
  /* Compute the difference of new - old */
  Grid delta = q_new - q_old;

  for (auto& P : Ps) {
    /* This is the PIC/FLIP update */
    P [ID] = alpha * q_new.lerp (P.p) + (1. - alpha) * (P [ID] + delta.lerp (P.p));
  }
}

void Particles::adjust (const VectorGrid& v, const Grid& rho, const LevelSet& fls, const LevelSet& sls) {
  recount ();
  prune (fls);
  seed (fls, sls, v, rho);
}

void Particles::advect (const VectorGrid& v, const LevelSet& sls, const double dt) {
  std::vector<size_t> dirty;

  size_t id = 0;
  for (auto& P : Ps) {
    P.p = RK3 (P.p, dt, v);
    Vec new_p = sls.find_closest (P.p);

    /* Check if good solution was found */
    if (new_p.x > 0 and new_p.y > 0) P.p = new_p;
    /* Otherwise simply discard particle */
    else dirty.push_back (id);

    id++;
  }

  for (const auto& id : dirty) {
    Ps [id] = Ps.back ();
    Ps.pop_back ();
  }
}

void Particles::recount () const {
  count.clear ();

  for (const auto& P : Ps) {
    int i = (int) P.p.x;
    int j = (int) P.p.y;

    count [id (i, j)]++;
  }
}

void Particles::prune (const LevelSet& fls) {
  /* Here we prune particles landing outside fluid or overcrowding area */
  for (auto it = Ps.begin (); it != Ps.end ();) {
    /* Check if position is outside liquid or too crowded */
    id idx ((int) (*it).p.x, (int) (*it).p.y);

    if (fls.lerp ((*it).p) > 0 or count [idx] > maxPerCell) {it = Ps.erase (it); count [idx] -= 1;}
    else ++it;
  }
}

void Particles::seed (const LevelSet& fls, const LevelSet& sls, const VectorGrid& v, const Grid& rho) {
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      int idx = i + j * w;

      /* Check if a new particles is needed */
      if (fls [idx] > 0) continue;

      for (int n = 0; n < minPerCell - count [idx]; n++) {
        /* Create new position for particles */
        double x = i + R.next ();
        double y = j + R.next ();

        Vec new_pos (x, y);

        /* Reject point if in solid */
        if (sls.lerp (new_pos) < 0.) continue;

        /* Reject point if is too far from fluid */
        if (fls.lerp (new_pos) > 0.25 * v.dx) continue;

        /* Adjust for small out of fluid position */
        if (fls.lerp (new_pos) > 0.) {new_pos = fls.find_closest (new_pos); assert (new_pos.x > 0);}

        Vec new_vel (v.x.lerp (new_pos), v.y.lerp (new_pos));
        std::array<double, 1> new_rho = {rho.lerp (new_pos)};

        Ps.push_back (Particle (new_pos, new_vel, new_rho));
      }
    }
  }
}
