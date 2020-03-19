#include <fstream>
#include <sstream>
#include <iomanip>

#include "include/solver.h"
#include "include/particles.h"
#include "include/lodepng.h"

int main (int argc, char** argv) {
  const int L = 1;
  const int w = 512;
  const int h = 512;
  const double dx = (double) (L) / w;
  const double dy = (double) (L) / h;
  const double dt = 0.1;
  const double r = 0.7;
  const double alpha = .01;

  unsigned char *img = new unsigned char [w * h * 4];

  Grid rho (w, h);
  LevelSet fls (w, h);
  LevelSet sls ("../res/solid_init.txt");

  VectorGrid v (w, h, dx, dy); //, "../res/vx_init.txt", "../res/vy_init.txt");
  VectorGrid g (w, h, dx, dy); //, "../res/fx_init.txt", "../res/fy_init.txt");

  /* Here we add the gravity force */
  g.y += 0.5;

  /* This is the pressure solver */
  Solver solver (w, h, dx, dy, dt, fls);

  // std::ofstream rho_out ("../out/rho.txt");
  // std::ofstream fls_out ("../out/fls.txt");
  // std::ofstream sls_out ("../out/sls.txt");
  // std::ofstream vx_out ("../out/vx.txt");
  // std::ofstream vy_out ("../out/vy.txt");
  // std::ofstream part_out ("../out/part.txt");

  Grid count (w, h);
  Particles P (r, count, "../res/fluid_init.txt");

  int T = argc > 1 ? std::stoi (argv [1]) : 50;

  for (int t = 0; t < T; t++) {
    std::cout << "Building Frame: " << t + 1 << "/" << T
              << " (Particles: " << P.Ps.size() << ")" << '\r' << std::flush;
    /* First we grab the particles values into the grid */
    fls.fit (P);
    fls &= sls;

    v.fit (P);
    rho.fit (P, 0);

    /* Here we prune and seed particles */
    // P.adjust (v, rho, fls, sls);

    /* Here we add the external forces */
    v += g * dt;


    /* Here we build the solid levelset */
    // For the moment is constant, we no not need to modify it.

    /* Then we solve for pressure */
    VectorGrid v_new = solver.project (fls, sls, v);

    P.fit (v, v_new, alpha);


    /* Finally we update the particles through the divergence free v-field */
    P.advect (v_new, sls, dt);

    /* Here we save current frame */
    std::ostringstream spath;
    spath << "../frames/" << std::setw (3) << std::setfill('0') << t << ".png" << std::setfill (' ');

    const std::string& path = spath.str ();

    img << fls.volume ();
    lodepng_encode32_file (path.c_str (), img, w, h);
  }

  std::cout << "\033[K\rSimulation Completed." << '\n';

  return 0;
}
