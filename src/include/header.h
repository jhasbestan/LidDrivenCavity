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
#include <unordered_map>
#include <vector>
#include <stdexcept>
using namespace std;

#define RESTRICT 1

//vector<bitset<M>> mesh;



class Q
{
private:


#if(1)

double * __restrict__ u=nullptr;
double * __restrict__ v=nullptr;
double * __restrict__ p=nullptr;

double * __restrict__ up=nullptr;
double * __restrict__ vp=nullptr;
double * __restrict__ pp=nullptr;

double * __restrict__ un=nullptr;
double * __restrict__ vn=nullptr;
double * __restrict__ pn=nullptr;

double * __restrict__ Res=nullptr;

#else

double * u=nullptr;
double *v=nullptr;
double *p=nullptr;

double *up=nullptr;
double *vp=nullptr;
double *pp=nullptr;

double *un=nullptr;
double *vn=nullptr;
double *pn=nullptr;

double * __restrict__ Res=nullptr;



#endif



double dx;
double dy=dx,dz=dx;
double dt;


uint sizeS;
uint sizeP;

uint longEnd;
uint shortEnd;

public:
Q(int nmax1,double dx,double dy);

int nmax;
int index(int i,int j,int k);

void initialize(int N,double *a);

void VTK_out( double *X, double *Y, int N );

void VTK_out_with_ghost( double *X, double *Y);

void Struct_2D_Ghost(double Xa,double Xb, double Ya, double Yb, double **X, double **Y);

void KovFlow(double x,double y,double *ux, double *uy,double *p1);

void setExactBC(double Xa,double Xb,double Ya,double Yb);

void Uxy(double Xa,double Ya,int i,int j,double *xy);
void Vxy(double Xa,double Ya,int i,int j,double *xy);
void Pxy(double Xa,double Ya,int i,int j,double *xy);
void Grad();
inline double U(int i,int j);

void setBoundaryLidDrivenCavity();

void showExact(double Xa,double Ya);
int uIdx( int i,int j);
int vIdx(int i,int j);
int wIdx(int i,int j);
int pIdx(int i,int j);

// i=0, face 0; i=1, face 1;
// j=0, face 2; i=1, face 3;
// i=0, face 4; i=1, face 5;
void uGetFace(uint fId,uint *size , double **face);
void vGetFace(uint fId,uint *size , double **face);
void wGetFace(uint fId,uint *size , double **face);

void getRes(double Re);
void getResTotal(double Re);
void predict();
void project();
void correct();
void update();
void start();
void debug(double Xa,double Ya);

void getResNorm(double *del_u);
void uSetFace( double *val );
void setNeumanPressure();
void setNeumanPressureLDC();
void setSolid(  );

~Q();
};





#endif
