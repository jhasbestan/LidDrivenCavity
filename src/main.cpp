#include "header.h"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#define LID 0

double rhs(double &x, double &y, double &z);
void Struct_2D(double Xa, double Xb, double Ya, double Yb, int64_t N, int64_t M,
               double **X, double **Y);
void VTK_out(double *X, double *Y, double *Z, int64_t N, int64_t M, int64_t ID);

int main(int argsc, char *argsv[]) {

  int64_t N=atoi(argsv[1]);
  int64_t count_max=atoi(argsv[2]);
  if(argsc<3)
  {
    cout << "forgot to pass in N and/or number of iter and/or number of iterss" << endl;
  }
  // N is the number if elements

  // include the ghost inside initilization
  double a[4] = {0.0, 0.0, 0.0, 0.0};

  double *d = nullptr;
/*
  double *X = nullptr;
  double *Y = nullptr;
*/
#ifndef LID
  double Xa = -.5, Xb = 1.0;
  double Ya = -.5, Yb = 0.5;
#else
  double Xa = 0.0, Xb = 1.0;
  double Ya = 0.0, Yb = 1.0;
#endif

  double dx = (Xb - Xa) / N;
  double dy = (Yb - Ya) / N;

  Q q(N, dx, dy);

  q.initialize(N, a);
  delete[] d;

  double Re = 40.0;

  int64_t count = 0;
  double res = 1.0;
  double total_time=0.0;

#ifndef LID
  q.setExactBC(Xa, Xb, Ya, Yb);
  q.start();
  count = 0;
  while (res > 1.e-7) {
    double start_time = omp_get_wtime();
    q.getRes(Re);
    q.predict();
    q.project();
    q.setNeumanPressure();
    q.correct();
    q.update();

    q.getResTotal(Re);
    q.getResNorm(&res);
    double end_time = omp_get_wtime();
    count++;
    if (count % 60 == 0) {
      cout << res << endl;
    }
    total_time+=end_time-start_time;
  }
//   std::cout<< " total time Kov "<< total_time <<std::endl;
#else
  // solve lid driven Cavity
  Re = 100.;

  q.setBoundaryLidDrivenCavity();
  q.start();
  count = 0;
  double res_old = 2.0;
#ifdef _OPENMP
  double start_time = omp_get_wtime();
#endif
//  while (res > 1.e-7) {
  while (count < count_max) {
    q.getRes(Re);
    q.predict();
    q.setNeumanPressureLDC();
    q.project();
    q.setBoundaryLidDrivenCavity();
    q.correct();
    q.update();
    res_old = res;
    q.getResTotal(Re);
    q.getResNorm(&res);
    count++;
/*
    if (count % 100 == 0) {
      cout << res << endl;
    }
*/
    if (fabs(res_old - res) < 1e-14) {
      // break;
    }
    res_old = res;
  }
#endif
#ifdef _OPENMP
    double end_time = omp_get_wtime();
    total_time += end_time-start_time;
#endif
    cout <<" residual = "<< res << endl;

/*
#if (1)
  Struct_2D(Xa, Xb, Ya, Yb, N, N, &X, &Y);
  q.VTK_out(X, Y, N);

#else
  q.Struct_2D_Ghost(Xa, Xb, Ya, Yb, &X, &Y);
  q.VTK_out_with_ghost(X, Y);
#endif

   std::cout<< " total time "<< total_time <<std::endl;
*/
};

void Struct_2D(double Xa, double Xb, double Ya, double Yb, int64_t N, int64_t M,
               double **X, double **Y) {

  double hx = N;
  double hy = N;
  double Xh = (Xb - Xa) / (hx);
  double Yh = (Yb - Ya) / (hy);

  (*X) = new double[N + 1];
  (*Y) = new double[N + 1];

  for (int64_t i = 0; i < N + 1; i++) {
    (*X)[i] = Xa + Xh * i;
  }

  for (int64_t i = 0; i < N + 1; i++) {
    (*Y)[i] = Ya + Yh * i;
  }
}

double rhs(double &x, double &y, double &z) {
  static const double pi = atan(1.0) * 4.0;
  static const double pi_2 = pi * pi;
  return (-(3.0 * pi_2 * sin((pi * x) / 2.0) * sin((pi * y) / 2.0) *
            sin((pi * z) / 2.0)) /
          4.0);
}
