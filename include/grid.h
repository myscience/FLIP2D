#ifndef GRID_H
#define GRID_H

#include "utils.h"

class VectorGrid;
class Particles;
class LevelSet;

class Grid {
public:
  Grid (const int& _size);
  Grid (const int& w, const int& h);
  Grid (const int& w, const int& h, const int offX, const int offY);

  Grid (std::string filename);
  Grid (const double offX, const double offY, std::string filename);
  ~Grid ();

  Grid (const Grid& a);

  double& operator[] (int idx);
  const double& operator[] (int idx) const;

  uint& operator() (int i);
  const uint& operator() (int i) const;

  double& operator[] (id idx);
  const  double& operator[] (id idx) const;

  uint& operator() (id idx);
  const uint& operator() (id idx) const;

  /* Assignment operators */
  Grid& operator= (const Grid& b);

  double operator*(const Grid& a) const;
  Grid operator+ (const Grid& b) const;
  Grid operator- (const Grid& b) const;
  Grid operator* (const double c);
  Grid operator- (const double b) const;
  Grid operator* (const double c) const;

  Grid& operator*= (const double c);
  Grid& operator+= (const double c);
  Grid& operator+= (const Grid& b);
  Grid& operator-= (const double c);
  Grid& operator-= (const Grid& b);

  bool operator== (double c) const;

  Grid smooth (int sx, int sy);
  double cerp (const Vec& pos) const;
  double lerp (const Vec& pos) const;
  Vec lerpGrad (const Vec& pos) const;

  double norm_inf () const;
  double norm_2 () const;

  void clear ();

  void advect (const double dt, const VectorGrid& v);
  void extrapolate ();
  void fit (const Particles& P, size_t ID);

  double _cerp (double a, double b, double c, double d, double x) const;
  double _lerp (double a, double b, double x) const;
  double _h (double r) const;
  double Bspline (const Vec& p, const id& idx);


  friend std::ostream& operator<< (std::ostream& stream, const Grid& grid);
  friend std::istream& operator>> (std::ifstream& stream, Grid& grid);

  friend unsigned char* operator<< (unsigned char* img, const Grid& grid);

  int w;
  int h;
  int size;

  double offX;
  double offY;

  double* data = nullptr;
  uint*  _data = nullptr;

  double nullvalue = 0.;
  uint nulllabel = INVALID;

private:
};

#endif
