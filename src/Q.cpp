#include "header.h"
#include <cmath>
#include <cstring>
#define ZSHRINK 0.1
#define AVG 1
// define poisson with FVM
#define PFV 0

Q::Q(int nmax1, double dx1, double dy1) {
  // nmax1 is the number of elements

  dx = dx1;
  dy = dy1;
  dt = 0.00001*dx;

  nmax = nmax1 + 2;

  cout << dx << " " << dy << " " << dz << " " << dt << endl;
  //
  longEnd = nmax;

  shortEnd = nmax - 1;

  // nmax is the number of elements with ghost cell

  cout << "nelem " << nmax1 << endl;

  cout << "nelem with ghost" << nmax << endl;

  sizeS = (nmax + 1) * (nmax);

  std::cout<< " Number of Elements " << sizeS <<std::endl;
  // u is u at current time step, un means new values hence at time step n+1 and
  // up is previous value of u means u at n-1
  /*
  u = new double[sizeS];
  v = new double[sizeS];

  up = new double[sizeS];
  vp = new double[sizeS];

  un = new double[sizeS];
  vn = new double[sizeS];
*/
  Res = new double[2 * sizeS];

  sizeP = nmax * nmax;

  posix_memalign((void**)&u,64,sizeS*sizeof(double));
  posix_memalign((void**)&up,64,sizeS*sizeof(double));
  posix_memalign((void**)&un,64,sizeS*sizeof(double));
  posix_memalign((void**)&v,64,sizeS*sizeof(double));
  posix_memalign((void**)&vp,64,sizeS*sizeof(double));
  posix_memalign((void**)&vn,64,sizeS*sizeof(double));
#ifdef __ICC
  __assume_aligned(u,64);
  __assume_aligned(up,64);
  __assume_aligned(un,64);
  __assume_aligned(v,64);
  __assume_aligned(vp,64);
  __assume_aligned(vn,64);
#endif
 // p = new double[sizeP];
//  pn = new double[sizeP];
  posix_memalign((void**)&p,64,sizeP*sizeof(double));
  posix_memalign((void**)&pn,64,sizeP*sizeof(double));
  posix_memalign((void**)&pn_old,64,sizeP*sizeof(double));
  posix_memalign((void**)&pp,64,sizeP*sizeof(double));
#ifdef __ICC
  __assume_aligned(pn_old,64);
  __assume_aligned(pn,64);
#endif
//  pn_old = new double[sizeP];
//  pp = new double[sizeP];

  // initialize the current time step

  double C = 0.0;
#pragma omp parallel for simd
  for (uint i = 0; i < sizeS; i++) {
    u[i] = C;
    v[i] = C;

    un[i] = C;
    vn[i] = C;

    up[i] = C;
    vp[i] = C;
  }

#pragma omp parallel for simd
  for (uint i = 0; i < 2 * sizeS; i++) {
    Res[i] = 0.0;
  }

#pragma omp parallel for simd
  for (uint i = 0; i < sizeP; i++) {
    p[i] = 0.0;
    pn[i] = 0.0;
    pp[i] = 0.0;
    pn_old[i]=0.0;
  }
}

Q::~Q() {
  free(u);
  free(v);
//  delete[] p;
  free(up);
  free(vp);
//  delete[] pp;
  free(un);
  free(vn);
//  delete[] pn;
 free(pn);
 free(pn_old);
 free(p);
 free(pp);
}

void Q::initialize(int N, double *val) {}

void Q::VTK_out(double *X, double *Y, int N) {

  FILE *fp = NULL;
  int M = N;

  // here we get some data into variable data
  char filename[64];
  sprintf(filename, "out%d.vtk", 0);
  fp = fopen(filename, "w");

  N++;
  M++;

  // fp=fopen("out.vtk","w");
  fprintf(fp, "# vtk DataFile Version 2.0 \n");
  fprintf(fp, "Grid\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET STRUCTURED_GRID\n");
  fprintf(fp, "DIMENSIONS %d %d %d\n", N, M, 1);
  fprintf(fp, "POINTS %d float\n", M * N);

    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        fprintf(fp, "%lf %lf %lf\n", X[i], Y[j], 0.0);
      }
  }

  fprintf(fp, "CELL_DATA %d\n", (N - 1) * (N - 1));

  fprintf(fp, "SCALARS U float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (int j = 1; j < N; j++) {
    for (int i = 1; i < N; i++) {
      fprintf(fp, "%lf\n", (u[uIdx(i, j)]));
    }
  }
  fprintf(fp, "SCALARS V float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (int j = 1; j < N; j++) {
    for (int i = 1; i < N; i++) {

      fprintf(fp, "%lf\n", (v[vIdx(i, j)]));
    }
  }
  fprintf(fp, "SCALARS P float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (int j = 1; j < N; j++) {
    for (int i = 1; i < N; i++) {
      fprintf(fp, " %lf\n", p[pIdx(i, j)]);
    }
  }
  fclose(fp);
}
#if (1)
void Q::VTK_out_with_ghost(double *X, double *Y) {
  FILE *fp = NULL;
  int M, N;

  // here we get some data into variable data
  char filename[64];
  sprintf(filename, "out%d.vtk", 0);
  fp = fopen(filename, "w");

  N = longEnd + 1;
  M = N;

  // fp=fopen("out.vtk","w");
  fprintf(fp, "# vtk DataFile Version 2.0 \n");
  fprintf(fp, "Grid\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET STRUCTURED_GRID\n");
  fprintf(fp, "DIMENSIONS %d %d %d\n", N, M, 1);
  fprintf(fp, "POINTS %d float\n", M * N);

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      fprintf(fp, "%lf %lf %lf\n", X[i], Y[j], 0.0);
    }
  }

  fprintf(fp, "CELL_DATA %d\n", (N - 1) * (N - 1));

  fprintf(fp, "SCALARS U float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (uint j = 0; j < longEnd; j++) {
    for (uint i = 0; i < longEnd; i++) {
#if (AVG)
      fprintf(fp, "%lf\n", (u[uIdx(i, j)] + u[uIdx(i + 1, j)]) * 0.5);
#else

      fprintf(fp, "%lf\n", (u[uIdx(i, j)]));
#endif
    }
  }

  fprintf(fp, "SCALARS V float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (unsigned int j = 0; j < longEnd; j++) {
    for (unsigned int i = 0; i < longEnd; i++) {
#if (AVG)
      fprintf(fp, "%lf\n", (v[vIdx(i, j)] + v[vIdx(i, j + 1)]) * 0.5);
#else
      fprintf(fp, "%lf\n", (v[vIdx(i, j)]));
#endif
    }
  }
  fprintf(fp, "SCALARS P float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");

  for (unsigned int j = 0; j < longEnd; j++) {
    for (unsigned int i = 0; i < longEnd; i++) {
      fprintf(fp, " %lf\n", p[pIdx(i, j)]);
    }
  }

  fclose(fp);
}
#endif

// indexes for different vector

inline int Q::uIdx(int i, int j) { return ((j) * (nmax + 1) + i); }

// i=0:nmax, j=0:nmax+1, k=0:nmax
inline int Q::vIdx(int i, int j) { return (j * (nmax) + i); }

// i=0:nmax, j=0:nmax, k=0:nmax+1
inline int Q::pIdx(int i, int j) { return (j * (nmax) + i); }

/************************************************
 *
 *
 ************************************************/

void Q::getRes(double Re) {

  double c1 = 0.5;
  double c2 = 0.5;

#pragma omp parallel 
{  
#pragma omp for nowait
    for (uint j = 1; j < shortEnd; j++) {
#pragma omp simd
  for (uint i = 1; i < longEnd; i++) {
      Res[uIdx(i, j)] =
          c1 * ((pow((u[uIdx(i + 1, j)] + u[uIdx(i, j)]), 2.) -
                 pow(u[uIdx(i, j)] + u[uIdx(i - 1, j)], 2.)) /
                    dx * 0.25 +
                0.25 / dy *
                    ((v[vIdx(i, j + 1)] + v[vIdx(i - 1, j + 1)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j + 1)]) -
                     (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j - 1)])))

          - 2. / Re / dx *
                ((u[uIdx(i + 1, j)] - u[uIdx(i, j)]) / dx -
                 (u[uIdx(i, j)] - u[uIdx(i - 1, j)]) / dx) -
          1. / Re / dy *
              (((u[uIdx(i, j + 1)] - u[uIdx(i, j)]) / dy +
                (v[vIdx(i, j + 1)] - v[vIdx(i - 1, j + 1)]) / dx) -
               ((u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx))

          // skew symmetric formulation
          // first trial dumb as sack of rocks, making analogy with Gaussian
          // Quadrature, need to use more points than one to integrate
          //
          /*
                           +c2*(u[uIdx(i,j)]*0.5/dx*(u[uIdx(i+1,j)]-u[uIdx(i,j)]+u[uIdx(i,j)]-u[uIdx(i-1,j)])
                           +0.25*(v[vIdx(i,j)]+v[vIdx(i,j+1)]+
             v[vIdx(i-1,j)]+v[vIdx(i-1,j+1)])*(u[uIdx(i,j+1)]-u[uIdx(i,j-1)])/2./dy);
          */
          + c2 * (0.25 * (u[uIdx(i, j)] + u[uIdx(i + 1, j)]) *
                      (u[uIdx(i + 1, j)] - u[uIdx(i, j)]) / dx +
                  0.25 * (u[uIdx(i, j)] + u[uIdx(i - 1, j)]) *
                      (u[uIdx(i, j)] - u[uIdx(i - 1, j)]) / dx +
                  (0.25 * (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                       (u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                   0.25 * (v[vIdx(i - 1, j + 1)] + v[vIdx(i, j + 1)]) *
                       (u[uIdx(i, j + 1)] - u[uIdx(i, j)]) / dy));

      // pressure grad
      //                  +1./dx*(p[pIdx(i,j)]-p[pIdx(i-1,j)]);
      //                     cout<<un[uIdx(i,j,k)]<<endl;
    }
  }
#pragma omp for 
    for (uint j = 1; j < longEnd; j++) {
#pragma omp simd
  for (uint i = 1; i < shortEnd; i++) {

      Res[sizeS + vIdx(i, j)] =
          c1 * ((pow((v[vIdx(i, j + 1)] + v[vIdx(i, j)]), 2.) -
                 pow(v[vIdx(i, j)] + v[vIdx(i, j - 1)], 2.)) /
                    dy * 0.25 +
                0.25 / dx *
                    ((v[vIdx(i + 1, j)] + v[vIdx(i, j)]) *
                         (u[uIdx(i + 1, j)] + u[uIdx(i + 1, j - 1)]) -
                     (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j - 1)]))) -
          2. / Re / dy *
              ((v[vIdx(i, j + 1)] - v[vIdx(i, j)]) / dy -
               (v[vIdx(i, j)] - v[vIdx(i, j - 1)]) / dy) -
          1. / Re / dx *
              (((u[uIdx(i + 1, j)] - u[uIdx(i + 1, j - 1)]) / dy +
                (v[vIdx(i + 1, j)] - v[vIdx(i, j)]) / dx) -
               ((u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx))

          // skew symmetric formulation
          /* dumb version
                           +c2*(v[vIdx(i,j)]*0.5/dy*(v[vIdx(i,j+1)]-v[vIdx(i,j)]+v[vIdx(i,j)]-v[vIdx(i,j-1)])
                           +0.25*(u[uIdx(i,j)]+u[uIdx(i+1,j)]+
             u[uIdx(i,j-1)]+u[uIdx(i+1,j-1)])*(v[vIdx(i+1,j)]-v[uIdx(i-1,j)])/2./dx);
          */

          + c2 * (0.25 * (v[vIdx(i, j)] + v[vIdx(i, j + 1)]) *
                      (v[vIdx(i, j + 1)] - v[vIdx(i, j)]) / dy +
                  0.25 * (v[vIdx(i, j)] + v[vIdx(i, j - 1)]) *
                      (v[vIdx(i, j)] - v[vIdx(i, j - 1)]) / dy +
                  (0.25 * (u[uIdx(i, j)] + u[uIdx(i, j - 1)]) *
                       (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx +
                   0.25 * (u[uIdx(i + 1, j - 1)] + u[uIdx(i + 1, j)]) *
                       (v[vIdx(i + 1, j)] - v[vIdx(i, j)]) / dx));

      // pressure grad
      //                 +1./dy*(p[pIdx(i,j)]-p[pIdx(i,j-1)]);

      //       cout<<vn[vIdx(i,j,k)]<<endl;
    }
  }
  }

  // cout<<"completed"<<endl;
}

void Q::getResTotal(double Re) {

  double c1 = 0.5;
  double c2 = 0.5;

#pragma omp parallel 
  {  
#pragma omp for nowait
    for (uint j = 1; j < shortEnd; j++) {
#pragma omp simd
  for (uint i = 1; i < longEnd; i++) {
      Res[uIdx(i, j)] =
          c1 * ((pow((u[uIdx(i + 1, j)] + u[uIdx(i, j)]), 2.) -
                 pow(u[uIdx(i, j)] + u[uIdx(i - 1, j)], 2.)) /
                    dx * 0.25 +
                0.25 / dy *
                    ((v[vIdx(i, j + 1)] + v[vIdx(i - 1, j + 1)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j + 1)]) -
                     (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j - 1)]))) -
          2. / Re / dx *
              ((u[uIdx(i + 1, j)] - u[uIdx(i, j)]) / dx -
               (u[uIdx(i, j)] - u[uIdx(i - 1, j)]) / dx) -
          1. / Re / dy *
              (((u[uIdx(i, j + 1)] - u[uIdx(i, j)]) / dy +
                (v[vIdx(i, j + 1)] - v[vIdx(i - 1, j + 1)]) / dx) -
               ((u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx))

          + c2 * (0.25 * (u[uIdx(i, j)] + u[uIdx(i + 1, j)]) *
                      (u[uIdx(i + 1, j)] - u[uIdx(i, j)]) / dx +
                  0.25 * (u[uIdx(i, j)] + u[uIdx(i - 1, j)]) *
                      (u[uIdx(i, j)] - u[uIdx(i - 1, j)]) / dx +
                  (0.25 * (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                       (u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                   0.25 * (v[vIdx(i - 1, j + 1)] + v[vIdx(i, j + 1)]) *
                       (u[uIdx(i, j + 1)] - u[uIdx(i, j)]) / dy));

      // pressure grad
      //+1. / dx *(p[pIdx(i, j)] - p[pIdx(i - 1, j)]);

      //                     cout<<un[uIdx(i,j,k)]<<endl;
    }
  }

#pragma omp for
    for (uint j = 1; j < longEnd; j++) {
#pragma omp simd
  for (uint i = 1; i < shortEnd; i++) {
      Res[sizeS + vIdx(i, j)] =
          c1 * ((pow((v[vIdx(i, j + 1)] + v[vIdx(i, j)]), 2.) -
                 pow(v[vIdx(i, j)] + v[vIdx(i, j - 1)], 2.)) /
                    dy * 0.25 +
                0.25 / dx *
                    ((v[vIdx(i + 1, j)] + v[vIdx(i, j)]) *
                         (u[uIdx(i + 1, j)] + u[uIdx(i + 1, j - 1)]) -
                     (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j - 1)]))) -
          2. / Re / dy *
              ((v[vIdx(i, j + 1)] - v[vIdx(i, j)]) / dy -
               (v[vIdx(i, j)] - v[vIdx(i, j - 1)]) / dy) -
          1. / Re / dx *
              (((u[uIdx(i + 1, j)] - u[uIdx(i + 1, j - 1)]) / dy +
                (v[vIdx(i + 1, j)] - v[vIdx(i, j)]) / dx) -
               ((u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx))

          + c2 * (0.25 * (v[vIdx(i, j)] + v[vIdx(i, j + 1)]) *
                      (v[vIdx(i, j + 1)] - v[vIdx(i, j)]) / dy +
                  0.25 * (v[vIdx(i, j)] + v[vIdx(i, j - 1)]) *
                      (v[vIdx(i, j)] - v[vIdx(i, j - 1)]) / dy +
                  (0.25 * (u[uIdx(i, j)] + u[uIdx(i, j - 1)]) *
                       (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx +
                   0.25 * (u[uIdx(i + 1, j - 1)] + u[uIdx(i + 1, j)]) *
                       (v[vIdx(i + 1, j)] - v[vIdx(i, j)]) / dx));

      // pressure grad
      //+1. / dy *(p[pIdx(i, j)] - p[pIdx(i, j - 1)]);
    }
  }
  }
}

void Q::start() {

  sizeS = (nmax + 1) * (nmax);

#pragma omp parallel for simd
  for (uint i = 0; i < sizeS; i++) {
    un[i] = u[i];
    up[i] = u[i];

    vp[i] = v[i];
    vn[i] = v[i];
  }
#pragma omp parallel for simd
  for (uint i = 0; i < sizeP; i++) {
    pn[i] = p[i];
    pp[i] = p[i];
  }
}

void Q::update() {

#pragma omp parallel for simd
  for (uint i = 0; i < sizeS; i++) {
    up[i] = u[i];
    vp[i] = v[i];
    u[i] = un[i];
    v[i] = vn[i];
  }

#pragma omp parallel for simd
  for (uint i = 0; i < sizeP; i++) {
    pp[i] = p[i];
    p[i] = pn[i];
  }
}

void Q::getResNorm(double *del_u) {

  double res0, res1;
  sizeS = (nmax + 1) * (nmax);
#pragma omp parallel for simd reduction(+:res0,res1)
  for (uint i = 0; i < sizeS; i++) {
    res0 +=  (un[i] - up[i]) * (un[i] - up[i]);
    res1 +=  (vn[i] - vp[i]) * (vn[i] - vp[i]);
  }
  *del_u = sqrt((res0 + res1) / 3. / sizeS);
}

#if (!PFV)
// FD formulation
void Q::project() {

  //std::cout<<" calling project ... "<<std::endl;
  // Gaus-Seidel Iteration
  double err=1.0;
  double val;

  double c1 = 2. / dx / dx + 2. / dy / dy;

  double *tmp;
  

#ifdef __ICC
__assume_aligned(pn_old,64);
__assume_aligned(pn,64);
#endif

//  pn[pIdx(5, 5)] = 0.0;
  //  while ( err > 1.e-12 )
  for (uint l = 0; l < 10; l++) {
#pragma omp parallel for 
      for (uint j = 1; j < shortEnd; j++) {
#pragma omp simd
    for (uint i = 1; i < shortEnd; i++) {
          val = -((un[uIdx(i + 1, j)] - un[uIdx(i, j)]) / dx +
                  (vn[vIdx(i, j + 1)] - vn[vIdx(i, j)]) / dy) /
                    dt +
                (pn_old[pIdx(i + 1, j)] + pn_old[pIdx(i - 1, j)]) / dx / dx +
                (pn_old[pIdx(i, j + 1)] + pn_old[pIdx(i, j - 1)]) / dy / dy;
          pn[pIdx(i, j)] = val / c1;

      }
    }
     tmp    = pn_old;
     pn_old = pn;
     pn     = tmp;
  }
}

#else

void Q::project() {
  std::cout<<" calling project "<<std::endl;
  // Gaus-Seidel Iteration
  double err = 1.0;
  double val;

  double c1 = 2. / dx / dx + 2. / dy / dy;

  //    while ( err > 1.e-12 )
  for (uint l = 0; l < 12; l++) {
    err = 0.0;
    pn[pIdx(0, 0)] = 0.0;
#pragma omp parallel for
    for (uint i = 1; i < shortEnd; i++) {
#pragma omp simd
      for (uint j = 1; j < shortEnd; j++) {
        val = ((un[uIdx(i + 1, j)] / dt -
                (pn[pIdx(i + 1, j)] - pn[pIdx(i, j)]) / dx) +
               -(un[uIdx(i, j)] / dt -
                 (pn[pIdx(i, j)] - pn[pIdx(i - 1, j)]) / dx)) /
                  dx +
              +((vn[vIdx(i, j + 1)] / dt -
                 (pn[pIdx(i, j + 1)] - pn[pIdx(i, j)]) / dy) +
                -(vn[vIdx(i, j)] / dt -
                  (pn[pIdx(i, j)] - pn[pIdx(i, j - 1)]) / dy)) /
                  dy;
    //    err += fabs(val);
        pn[pIdx(i, j)] = pn[pIdx(i, j)] - val / c1;
        //                               cout<<"val "<<val<<endl;
      }
    }
  }
}
#endif

void Q::predict() {

#pragma omp parallel for simd
  for (uint i = 0; i < sizeS; i++) {
    un[i] = -Res[i] * dt + u[i];
    vn[i] = -Res[sizeS + i] * dt + v[i];
  }
}

void Q::correct() {

#pragma omp parallel   
{
#pragma omp for nowait  
    for (uint j = 1; j < shortEnd; j++) {
#pragma omp simd
  for (uint i = 1; i < longEnd; i++) {
      un[uIdx(i, j)] =
          -dt * (pn[pIdx(i, j)] - pn[pIdx(i - 1, j)]) / dx + un[uIdx(i, j)];
    }
  }
#pragma omp for 
    for (uint j = 1; j < longEnd; j++) {
#pragma omp simd
  for (uint i = 1; i < shortEnd; i++) {
      vn[vIdx(i, j)] =
          -dt * (pn[pIdx(i, j)] - pn[pIdx(i, j - 1)]) / dy + vn[vIdx(i, j)];
    }
  }
}
}

// this is for solid

void Q::setNeumanPressure() {
  // apply newman bc on the solid
  // we solid faces so five neuman or solid
#pragma omp parallel for 
  for (int i = 0; i < nmax; i++) {
    pn[pIdx(i, nmax - 1)] = pn[pIdx(i, nmax - 2)];
    pn[pIdx(i, 0)] = pn[pIdx(i, 1)];
  }
  //        pn[(pIdx(1,1,1))]=0.0;
}

#if (0)
void Q::uSetFace(double *val) {
  for (uint j = 0; j < nmax + 1; j++) {
    for (uint k = 0; k < nmax; k++) {
      v[vIndex(0, j, k)] = (*val) * 2.0 - v[vIndex(1, j, k)];

      vp[vIndex(0, j, k)] = (*val) * 2.0 - v[vIndex(1, j, k)];

      vn[vIndex(0, j, k)] = (*val) * 2.0 - v[vIndex(1, j, k)];
    }
  }

  //
  //
  // set u and w as zero
  //
  //
  for (uint j = 0; j < nmax; j++) {
    for (uint k = 0; k < nmax; k++) {
      u[uIndex(0, j, k)] = 0;
      up[uIndex(0, j, k)] = 0;
      un[uIndex(0, j, k)] = 0;

      u[uIndex(nmax, j, k)] = 0;
      up[uIndex(nmax, j, k)] = 0;
      un[uIndex(nmax, j, k)] = 0;
    }
  }
  for (uint j = 0; j < nmax + 1; j++) {
    for (uint k = 0; k < nmax; k++) {
      u[uIndex(j, 0, k)] = 0;
      up[uIndex(j, 0, k)] = 0;
      un[uIndex(j, 0, k)] = 0;

      u[uIndex(j, nmax - 1, k)] = 0;
      up[uIndex(j, nmax - 1, k)] = 0;
      un[uIndex(j, nmax - 1, k)] = 0;

      u[uIndex(j, k, 0)] = 0;
      up[uIndex(j, k, 0)] = 0;
      un[uIndex(j, k, 0)] = 0;

      u[uIndex(j, k, nmax - 1)] = 0;
      up[uIndex(j, k, nmax - 1)] = 0;
      un[uIndex(j, k, nmax - 1)] = 0;
    }
  }

  //
  //
  //
  //

  for (uint j = 0; j < nmax; j++) {
    for (uint k = 0; k < nmax + 1; k++) {
      w[wIndex(0, j, k)] = 0;
      wp[wIndex(0, j, k)] = 0;
      wn[wIndex(0, j, k)] = 0;

      w[wIndex(nmax - 1, j, k)] = 0;
      wp[wIndex(nmax - 1, j, k)] = 0;
      wn[wIndex(nmax - 1, j, k)] = 0;

      w[wIndex(j, 0, k)] = 0;
      wp[wIndex(j, 0, k)] = 0;
      wn[wIndex(j, 0, k)] = 0;

      w[wIndex(j, nmax - 1, k)] = 0;
      wp[wIndex(j, nmax - 1, k)] = 0;
      wn[wIndex(j, nmax - 1, k)] = 0;
    }
  }

  for (uint j = 0; j < nmax; j++) {
    for (uint k = 0; k < nmax; k++) {
      w[wIndex(j, k, 0)] = 0;
      wp[wIndex(j, k, 0)] = 0;
      wn[wIndex(j, k, 0)] = 0;

      w[wIndex(j, k, nmax - 1)] = 0;
      wp[wIndex(j, k, nmax - 1)] = 0;
      wn[wIndex(j, k, nmax - 1)] = 0;
    }
  }
}

//
// enforce solid
//

void Q::setSolid() {

  //
  //
  // set u and w  and v negative value of the closest interior
  //
  // x=-1 and x=1. planes
  //
  //
  //
  //
  for (uint j = 0; j < nmax; j++) {
    for (uint k = 0; k < nmax; k++) {

      u[uIndex(1, j, k)] = 0.0;
      up[uIndex(1, j, k)] = 0.0;
      un[uIndex(1, j, k)] = 0.0;

      u[uIndex(0, j, k)] = -u[uIndex(2, j, k)];
      up[uIndex(0, j, k)] = -up[uIndex(2, j, k)];
      un[uIndex(0, j, k)] = -un[uIndex(2, j, k)];

      u[uIndex(nmax, j, k)] = -u[uIndex(nmax - 2, j, k)];
      up[uIndex(nmax, j, k)] = -up[uIndex(nmax - 2, j, k)];
      un[uIndex(nmax, j, k)] = -un[uIndex(nmax - 2, j, k)];

      u[uIndex(nmax - 1, j, k)] = 0.0;
      up[uIndex(nmax - 1, j, k)] = 0.0;
      un[uIndex(nmax - 1, j, k)] = 0.0;
    }
  }

  for (uint j = 1; j < nmax; j++) {
    for (uint k = 1; k < nmax - 1; k++) {
      /*
            v[vIndex( 0, j, k )] = - v[vIndex( 1, j, k )];
            vp[vIndex( 0, j, k )] = - vp[vIndex( 1, j, k )];
            vn[vIndex( 0, j, k )] = - vn[vIndex( 1, j, k )] ;
     */
      v[vIndex(nmax - 1, j, k)] = -v[vIndex(nmax - 2, j, k)];
      vp[vIndex(nmax - 1, j, k)] = -vp[vIndex(nmax - 2, j, k)];
      vn[vIndex(nmax - 1, j, k)] = -vn[vIndex(nmax - 2, j, k)];
    }
  }

  for (uint j = 1; j < nmax - 1; j++) {
    for (uint k = 1; k < nmax; k++) {

      w[wIndex(0, j, k)] = -w[wIndex(1, j, k)];
      wp[wIndex(0, j, k)] = -wp[wIndex(1, j, k)];
      wn[wIndex(0, j, k)] = -wn[wIndex(1, j, k)];

      w[wIndex(nmax - 1, j, k)] = -w[wIndex(nmax - 2, j, k)];
      wp[wIndex(nmax - 1, j, k)] = -wp[wIndex(nmax - 2, j, k)];
      wn[wIndex(nmax - 1, j, k)] = -wn[wIndex(nmax - 2, j, k)];
    }
  }

  // ==============================================================
  //
  //   y=-1, y=+1 faces
  //
  // =============================================================
  //

  for (uint j = 1; j < nmax; j++) {
    for (uint k = 1; k < nmax - 1; k++) {

      u[uIndex(j, 0, k)] = -u[uIndex(j, 1, k)];
      up[uIndex(j, 0, k)] = -up[uIndex(j, 1, k)];
      un[uIndex(j, 0, k)] = -un[uIndex(j, 1, k)];

      u[uIndex(j, nmax - 1, k)] = -u[uIndex(nmax - 2, j, k)];
      up[uIndex(j, nmax - 1, k)] = -up[uIndex(nmax - 2, j, k)];
      un[uIndex(j, nmax - 1, k)] = -un[uIndex(nmax - 2, j, k)];
    }
  }

  for (uint j = 1; j < nmax; j++) {
    for (uint k = 1; k < nmax; k++) {

      v[vIndex(j, nmax, k)] = -v[vIndex(j, nmax - 1, k)];
      vp[vIndex(j, nmax, k)] = -vp[vIndex(j, nmax - 1, k)];
      vn[vIndex(j, nmax, k)] = -vn[vIndex(j, nmax - 1, k)];

      v[vIndex(j, 1, k)] = 0.0;
      vp[vIndex(j, 1, k)] = 0.0;
      vn[vIndex(j, 1, k)] = 0.0;

      v[vIndex(j, nmax - 1, k)] = 0.0;
      vp[vIndex(j, nmax - 1, k)] = 0.0;
      vn[vIndex(j, nmax - 1, k)] = 0.0;

      v[vIndex(j, 0, k)] = -v[vIndex(j, 2, k)];
      vp[vIndex(j, 0, k)] = -vp[vIndex(j, 2, k)];
      vn[vIndex(j, 0, k)] = -vn[vIndex(j, 2, k)];
    }
  }

  for (uint j = 1; j < nmax - 1; j++) {
    for (uint k = 1; k < nmax; k++) {

      w[wIndex(j, 0, k)] = -w[wIndex(j, 1, k)];
      wp[wIndex(j, 0, k)] = -wp[wIndex(j, 1, k)];
      wn[wIndex(j, 0, k)] = -wn[wIndex(j, 1, k)];

      w[wIndex(j, nmax - 1, k)] = -w[wIndex(j, nmax - 2, k)];
      wp[wIndex(j, nmax - 1, k)] = -wp[wIndex(j, nmax - 2, k)];
      wn[wIndex(j, nmax - 1, k)] = -wn[wIndex(j, nmax - 2, k)];
    }
  }

  // set periodic boundary condition at Z direction

  for (uint j = 0; j < nmax; j++) {
    for (uint k = 0; k < nmax; k++) {

      w[wIndex(j, k, nmax - 1)] = w[wIndex(j, k, 0)];
      wp[wIndex(j, k, nmax - 1)] = wp[wIndex(j, k, 0)];
      wn[wIndex(j, k, nmax - 1)] = wn[wIndex(j, k, 0)];

      v[vIndex(j, k, nmax - 1)] = v[vIndex(j, k, 0)];
      vp[vIndex(j, k, nmax - 1)] = vp[vIndex(j, k, 0)];
      vn[vIndex(j, k, nmax - 1)] = vn[vIndex(j, k, 0)];

      u[uIndex(j, k, nmax - 1)] = u[uIndex(j, k, 0)];
      up[uIndex(j, k, nmax - 1)] = up[uIndex(j, k, 0)];
      un[uIndex(j, k, nmax - 1)] = un[uIndex(j, k, 0)];
    }
  }

  // periodic for pressure
  //
  for (uint j = 0; j < nmax; j++) {
    for (uint k = 0; k < nmax; k++) {

      pn[pIndex(j, k, nmax - 1)] = pn[pIndex(j, k, 0)];
    }
  }

  for (uint j = 1; j < nmax - 1; j++) {
    for (uint k = 1; k < nmax - 1; k++) {

      // y =-1 and +1 directions
      v[vIndex(j, nmax, k)] = -v[vIndex(j, nmax - 1, k)];
      vp[vIndex(j, nmax, k)] = -vp[vIndex(j, nmax - 1, k)];
      vn[vIndex(j, nmax, k)] = -vn[vIndex(j, nmax - 1, k)];

      v[vIndex(j, 0, k)] = -v[vIndex(j, 1, k)];
      vp[vIndex(j, 0, k)] = -vp[vIndex(j, 1, k)];
      vn[vIndex(j, 0, k)] = -vn[vIndex(j, 1, k)];

      // z=-1 and 1 directions

      w[wIndex(j, k, nmax)] = w[wIndex(j, k, 0)];
      wp[wIndex(j, k, nmax)] = wp[wIndex(j, k, 0)];
      wn[wIndex(j, k, nmax)] = wn[wIndex(j, k, 0)];
    }
  }

  for (uint j = 0; j < nmax + 1; j++) {
    for (uint k = 0; k < nmax; k++) {
      u[uIndex(j, 0, k)] = -u[uIndex(j, 1, k)];
      up[uIndex(j, 0, k)] = -up[uIndex(j, 1, k)];
      un[uIndex(j, 0, k)] = -un[uIndex(j, 1, k)];

      u[uIndex(j, nmax - 1, k)] = -u[uIndex(j, nmax - 2, k)];
      up[uIndex(j, nmax - 1, k)] = -up[uIndex(j, nmax - 2, k)];
      un[uIndex(j, nmax - 1, k)] = -un[uIndex(j, nmax - 2, k)];

      v[vIndex(nmax - 1, j, k)] = -v[vIndex(nmax - 2, j, k)];
      vp[vIndex(nmax - 1, j, k)] = -vp[vIndex(nmax - 2, j, k)];
      vn[vIndex(nmax - 1, j, k)] = -vn[vIndex(nmax - 2, j, k)];

      w[wIndex(0, k, j)] = -w[wIndex(1, k, j)];
      wp[wIndex(0, k, j)] = -wp[wIndex(1, k, j)];
      wn[wIndex(0, k, j)] = -wn[wIndex(1, k, j)];

      w[wIndex(nmax - 1, k, j)] = -w[wIndex(nmax - 2, k, j)];
      wp[wIndex(nmax - 1, k, j)] = -wp[wIndex(nmax - 2, k, j)];
      wn[wIndex(nmax - 1, k, j)] = -wn[wIndex(nmax - 2, k, j)];

      w[wIndex(k, 0, j)] = -w[wIndex(k, 1, j)];
      wp[wIndex(k, 0, j)] = -wp[wIndex(k, 1, j)];
      wn[wIndex(0, k, j)] = -wn[wIndex(k, 1, j)];

      w[wIndex(nmax - 1, k, j)] = -w[wIndex(nmax - 2, k, j)];
      wp[wIndex(nmax - 1, k, j)] = -wp[wIndex(nmax - 2, k, j)];
      wn[wIndex(nmax - 1, k, j)] = -wn[wIndex(nmax - 2, k, j)];
    }
  }
}

#endif
void Q::Struct_2D_Ghost(double Xa, double Xb, double Ya, double Yb, double **X,
                        double **Y)

{
  int i;

  double hx = nmax - 2;
  double hy = nmax - 2;
  double Xh = (Xb - Xa) / (hx);
  double Yh = (Yb - Ya) / (hy);

  cout << "delx " << Xh << endl;
  cout << "dely " << Yh << endl;

  (*X) = new double[nmax + 1];
  (*Y) = new double[nmax + 1];

  Xa = Xa - Xh;

  Ya = Ya - Yh;

#pragma omp parallel for simd 
  for (i = 0; i < nmax + 1; i++) {
    (*X)[i] = Xa + Xh * i;
  }

#pragma omp parallel for simd 
  for (i = 0; i < nmax + 1; i++) {
    (*Y)[i] = Ya + Yh * i;
  }
}

void Q::KovFlow(double x, double y, double *ux, double *uy, double *p1) {

  double Re = 20.0;
  double pi = atan(1.) * 4.;
  double lamda = Re / 2. - sqrt(Re * Re / 4.0 + 4. * pi * pi);

  *ux = 1. - exp(lamda * (x)) * cos(2. * pi * y);
  *uy = lamda / 2. / pi * exp(lamda * x) * sin(2. * pi * y);
  *p1 = 0.5 * (1.0 - exp(2. * lamda * x));
}

void Q::setExactBC(double Xa, double Xb, double Ya, double Yb) {
  double x1, y1, ux, uy, p1;
  double *xy = new double[2];

  for (int i = 0; i < nmax + 1; i++) {

    // set U on bottom and top
    Uxy(Xa, Ya, i, 0, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    u[uIdx(i, 0)] = ux;

    Uxy(Xa, Ya, i, nmax - 1, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    u[uIdx(i, nmax - 1)] = ux;

    // set v on sides
    // //
    //
    Vxy(Xa, Ya, 0, i, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    v[vIdx(0, i)] = uy;

    Vxy(Xa, Ya, nmax - 1, i, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    v[vIdx(nmax - 1, i)] = uy;
  }

  for (int i = 0; i < nmax; i++) {

    // set U on sides
    Uxy(Xa, Ya, 0, i, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    u[uIdx(0, i)] = ux;

    Uxy(Xa, Ya, nmax, i, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    u[uIdx(nmax, i)] = ux;

    // set v on top and bottom
    // //
    //
    Vxy(Xa, Ya, i, 0, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    v[vIdx(i, 0)] = uy;

    Vxy(Xa, Ya, i, nmax, xy);
    x1 = xy[0];
    y1 = xy[1];
    KovFlow(x1, y1, &ux, &uy, &p1);
    v[vIdx(i, nmax)] = uy;
  }
}

void Q::Uxy(double Xa, double Ya, int i, int j, double *xy) {

  double X1 = Xa - dx;
  xy[0] = X1 + i * dx;
  double Y1 = Ya - dy * 0.5;
  xy[1] = Y1 + j * dy;
}

void Q::Vxy(double Xa, double Ya, int i, int j, double *xy) {

  double X1 = Xa - dx * 0.5;
  xy[0] = X1 + i * dx;
  double Y1 = Ya - dy;
  xy[1] = Y1 + j * dy;
}

void Q::Pxy(double Xa, double Ya, int i, int j, double *xy) {

  double X1 = Xa - dx * 0.5;
  xy[0] = X1 + i * dx;
  double Y1 = Ya - 0.5 * dy;
  xy[1] = Y1 + j * dy;
}

void Q::showExact(double Xa, double Ya) {
  double x1, y1, ux, uy, p1;
  double *xy = new double[2];

  for (int i = 0; i < nmax + 1; i++) {
    for (int j = 0; j < nmax; j++) {

      Uxy(Xa, Ya, i, j, xy);
      x1 = xy[0];
      y1 = xy[1];
      KovFlow(x1, y1, &ux, &uy, &p1);
      u[uIdx(i, j)] = ux;
      // cout<<x1<<" "<<y1<<" "<<ux<<endl;
    }
  }
  for (int i = 0; i < nmax; i++) {
    for (int j = 0; j < nmax + 1; j++) {

      Vxy(Xa, Ya, i, j, xy);
      x1 = xy[0];
      y1 = xy[1];
      KovFlow(x1, y1, &ux, &uy, &p1);
      v[vIdx(i, j)] = uy;
      // cout<<x1<<" "<<y1<<" "<<ux<<endl;
    }
  }

  for (int i = 0; i < nmax; i++) {
    for (int j = 0; j < nmax; j++) {

      Pxy(Xa, Ya, i, j, xy);
      x1 = xy[0];
      y1 = xy[1];
      KovFlow(x1, y1, &ux, &uy, &p1);
      p[vIdx(i, j)] = p1;
      // cout<<x1<<" "<<y1<<" "<<ux<<endl;
    }
  }

  delete[] xy;
}

void Q::Grad() {

  for (int i = 0; i < nmax; i++) {
    for (int j = 0; j < nmax; j++) {

      // pn[pIdx(i,j)]=(un[uIdx(i+1,j)]-un[uIdx(i,j)])/dx;
      // pn[pIdx(i,j)]=(vn[vIdx(i,j+1)]-vn[vIdx(i,j)])/dy;
    }
  }

  for (int i = 0; i < nmax; i++) {
    for (int j = 0; j < nmax - 1; j++) {

      // pn[pIdx(i,j)]=(un[uIdx(i+1,j)]-un[uIdx(i,j)])/dx;
      // pn[pIdx(i,j)]=(un[uIdx(i,j+1)]-un[uIdx(i,j)])/dy;
    }
  }

  for (int i = 0; i < nmax - 1; i++) {
    for (int j = 0; j < nmax; j++) {

      // pn[pIdx(i,j)]=(un[uIdx(i+1,j)]-un[uIdx(i,j)])/dx;
      pn[pIdx(i, j)] = (vn[vIdx(i + 1, j)] - vn[vIdx(i, j)]) / dx;
    }
  }
}

double Q::U(int i, int j) { return (u[uIdx(i, j)]); }

void Q::debug(double Xa, double Ya) {

  int i = 4;
  int j = 5;

  double xy[2];

  Uxy(Xa, Ya, i, j, xy);
  cout << "U(i,j) coord " << xy[0] << " " << xy[1] << " value " << u[uIdx(i, j)]
       << endl;
  Uxy(Xa, Ya, i + 1, j, xy);
  cout << "U(i+1,j) coord " << xy[0] << " " << xy[1] << " value "
       << u[uIdx(i + 1, j)] << endl;
  Uxy(Xa, Ya, i, j + 1, xy);
  cout << "U(i,j+1) coord " << xy[0] << " " << xy[1] << " value "
       << u[uIdx(i, j + 1)] << endl;
  Uxy(Xa, Ya, i + 1, j + 1, xy);
  cout << "U(i+1,j+1) coord " << xy[0] << " " << xy[1] << " value "
       << u[uIdx(i + 1, j + 1)] << endl;
  Uxy(Xa, Ya, i - 1, j + 1, xy);
  cout << "U(i-1,j+1) coord " << xy[0] << " " << xy[1] << " value "
       << u[uIdx(i - 1, j + 1)] << endl;

  cout << "**********************************************************" << endl;
  Vxy(Xa, Ya, i, j, xy);
  cout << "V(i,j) coord " << xy[0] << " " << xy[1] << " value " << v[vIdx(i, j)]
       << endl;
  Vxy(Xa, Ya, i - 1, j, xy);
  cout << "V(i-1,j) coord " << xy[0] << " " << xy[1] << " value "
       << v[vIdx(i - 1, j)] << endl;
  Vxy(Xa, Ya, i, j + 1, xy);
  cout << "V(i,j+1) coord " << xy[0] << " " << xy[1] << " value "
       << v[vIdx(i, j + 1)] << endl;
  Vxy(Xa, Ya, i - 1, j + 1, xy);
  cout << "V(i-1,j+1) coord " << xy[0] << " " << xy[1] << " value "
       << v[vIdx(i - 1, j + 1)] << endl;
  Vxy(Xa, Ya, i - 1, j + 1, xy);
  cout << "V(i-1,j+1) coord " << xy[0] << " " << xy[1] << " value "
       << v[vIdx(i - 1, j + 1)] << endl;

  Vxy(Xa, Ya, i, j, xy);

  cout << "V coord " << xy[0] << " " << xy[1] << endl;

  Pxy(Xa, Ya, i, j, xy);

  double Re = 20.0;

  cout << "P coord " << xy[0] << " " << xy[1] << endl;

  double duu_dx, duv_dy, dp_dx;

  // duu_dx= -0.25*(pow( ( u[uIdx( i + 1, j )] + u[uIdx( i, j )] ), 2. ) - pow(
  // u[uIdx( i, j )] + u[uIdx( i - 1, j )], 2. ) ) / dx;
  duu_dx = pow((u[uIdx(i + 1, j)] + u[uIdx(i, j)]), 2.) -
           pow(u[uIdx(i, j)] + u[uIdx(i - 1, j)], 2.);

  duv_dy = 0.25 / dy *
           ((v[vIdx(i, j + 1)] + v[vIdx(i - 1, j + 1)]) *
                (u[uIdx(i, j)] + u[uIdx(i, j + 1)]) -
            (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                (u[uIdx(i, j)] + u[uIdx(i, j - 1)]));

  dp_dx = -1. / dx * (p[pIdx(i, j)] - p[pIdx(i - 1, j)]);

  cout << "duu_dx " << duu_dx << endl;

  cout << "duv_dy " << duv_dy << endl;

  cout << "dp_dx " << dp_dx << endl;

  double valx = -(pow((u[uIdx(i + 1, j)] + u[uIdx(i, j)]), 2.) -
                  pow(u[uIdx(i, j)] + u[uIdx(i - 1, j)], 2.)) /
                    dx * 0.25 -
                0.25 / dy *
                    ((v[vIdx(i, j + 1)] + v[vIdx(i - 1, j + 1)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j + 1)]) -
                     (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j - 1)])) -
                1. / dx * (p[pIdx(i, j)] - p[pIdx(i - 1, j)]) +
                2. / Re / dx *
                    ((u[uIdx(i + 1, j)] - u[uIdx(i, j)]) / dx -
                     (u[uIdx(i, j)] - u[uIdx(i - 1, j)]) / dx) +
                1. / Re / dy *
                    (((u[uIdx(i, j + 1)] - u[uIdx(i, j)]) / dy +
                      (v[vIdx(i, j + 1)] - v[vIdx(i - 1, j + 1)]) / dx) -
                     ((u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                      (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx));

  double valy = -(pow((v[vIdx(i, j + 1)] + v[vIdx(i, j)]), 2.) -
                  pow(v[vIdx(i, j)] + v[vIdx(i, j - 1)], 2.)) /
                    dy * 0.25 -
                0.25 / dx *
                    ((v[vIdx(i + 1, j)] + v[vIdx(i, j)]) *
                         (u[uIdx(i + 1, j)] + u[uIdx(i + 1, j - 1)]) -
                     (v[vIdx(i, j)] + v[vIdx(i - 1, j)]) *
                         (u[uIdx(i, j)] + u[uIdx(i, j - 1)])) -
                1. / dy * (p[pIdx(i, j)] - p[pIdx(i, j - 1)]) +
                2. / Re / dy *
                    ((v[vIdx(i, j + 1)] - v[vIdx(i, j)]) / dy -
                     (v[vIdx(i, j)] - v[vIdx(i, j - 1)]) / dy) +
                1. / Re / dx *
                    (((u[uIdx(i + 1, j)] - u[uIdx(i + 1, j - 1)]) / dy +
                      (v[vIdx(i + 1, j)] - v[vIdx(i, j)]) / dx) -
                     ((u[uIdx(i, j)] - u[uIdx(i, j - 1)]) / dy +
                      (v[vIdx(i, j)] - v[vIdx(i - 1, j)]) / dx));

  cout << "valx " << valx << endl;
  cout << "valy " << valy << endl;

  //   cout << "sizeS " << sizeS << endl;
}

void Q::setBoundaryLidDrivenCavity() {

#pragma omp parallel for simd
  for (int i = 1; i < nmax; i++) {
    // cond on U at i direction
    un[uIdx(i, 0)] = -un[uIdx(i, 1)];
    //        u[uIdx( i, nmax - 1 )] = 2. -  u[uIdx( i, nmax -2 )] ;
    un[uIdx(i, nmax - 1)] = 1.0;
    // cond on V at j direction
    vn[vIdx(0, i)] = -vn[vIdx(1, i)];
    vn[vIdx(nmax - 1, i)] = -vn[vIdx(nmax - 2, i)];
  }

  // no leaking at the top so the u velocity is zero
#pragma omp parallel for simd
  for (int i = 1; i < nmax - 1; i++) {
    // condition on V at i-direction
    //         vn[vIdx(i,1)]=0.0;
    //         vn[vIdx( i,0 )] = -vn[vIdx( 2, i )];
    // v[vIdx( i, nmax )] = 0.0;
    //         vn[vIdx(i,nmax)]=-vn[vIdx(i,nmax-1)];

    vn[vIdx(i, 1)] = 0.0;

    // condition on U in j direction
    //
    //        un[uIdx( 0, i )] = -un[uIdx( 2, i )];
    un[uIdx(1, i)] = 0.0;

    //        un[uIdx( nmax, i )] = -un[uIdx( nmax - 2, i )];
    un[uIdx(nmax - 1, i)] = 0.0;
  }
}

void Q::setNeumanPressureLDC() {
  // apply newman bc on the solid
  // we solid faces so five neuman or solid

#pragma omp parallel for simd
  for (int i = 0; i < nmax; i++) {
      pn[pIdx(i, 0)] = pn[pIdx(i, 1)];
      pn[pIdx(0, i)] = pn[pIdx(1, i)];
      pn[pIdx(nmax - 1, i)] = pn[pIdx(nmax - 2, i)];
  }
  //            pn[pIdx(0,0)]=0.0;
}
