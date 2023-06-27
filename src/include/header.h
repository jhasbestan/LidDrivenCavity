#ifndef _HEADER_H_
#define _HEADER_H_
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using namespace std;

#define RESTRICT 1




// vector<bitset<M>> mesh;

class Q {
private:

  double * __restrict__ u = nullptr;
  double * __restrict__ v = nullptr;
  double * __restrict__ p = nullptr;

  double * __restrict__ up = nullptr;
  double * __restrict__ vp = nullptr;
  double * __restrict__ pp = nullptr;

  double * __restrict__ un = nullptr;
  double * __restrict__ vn = nullptr;
  double * pn = nullptr;
  double * pn_old = nullptr;

  double * __restrict__ Res = nullptr;

  double dx;
  double dy;
  double dz;
  double dt;

  int64_t sizeS;
  int64_t sizeP;

  int64_t longEnd;
  int64_t shortEnd;

public:
  Q(int64_t nmax1, double dx, double dy);

  int64_t nmax;
  int64_t index(int64_t i, int64_t j, int64_t k);

  void initialize(int64_t N, double *a);

  void VTK_out(double *X, double *Y, int64_t N);

  void VTK_out_with_ghost(double *X, double *Y);

  void Struct_2D_Ghost(double Xa, double Xb, double Ya, double Yb, double **X,
                       double **Y);

  void KovFlow(double x, double y, double *ux, double *uy, double *p1);

  void setExactBC(double Xa, double Xb, double Ya, double Yb);

  void Uxy(double Xa, double Ya, int64_t i, int64_t j, double *xy);
  void Vxy(double Xa, double Ya, int64_t i, int64_t j, double *xy);
  void Pxy(double Xa, double Ya, int64_t i, int64_t j, double *xy);
  void Grad();
  inline double U(int64_t i, int64_t j);

  void setBoundaryLidDrivenCavity();

  void showExact(double Xa, double Ya);
  int64_t uIdx(int64_t i, int64_t j);
  int64_t vIdx(int64_t i, int64_t j);
  int64_t wIdx(int64_t i, int64_t j);
  int64_t pIdx(int64_t i, int64_t j);

  // i=0, face 0; i=1, face 1;
  // j=0, face 2; i=1, face 3;
  // i=0, face 4; i=1, face 5;
  void uGetFace(int64_t fId, int64_t *size, double **face);
  void vGetFace(int64_t fId, int64_t *size, double **face);
  void wGetFace(int64_t fId, int64_t *size, double **face);

  void getRes(double Re);
  void getResTotal(double Re);
  void predict();
  void project();
  void correct();
  void update();
  void start();
  void debug(double Xa, double Ya);

  void getResNorm(double *del_u);
  void uSetFace(double *val);
  void setNeumanPressure();
  void setNeumanPressureLDC();
  void setSolid();

  ~Q();
};

#endif
