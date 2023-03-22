#include "header.h"
#include <cmath>
#define LID 1

double rhs(double &x, double &y, double &z);
void Struct_2D(double Xa,double Xb, double Ya, double Yb,int N, int M,double **X, double **Y);
void VTK_out(double *X, double*Y,double *Z,int N, int M,int ID);

int main()
{

int N;
// N is the number if elements
cout<<"enter N"<<endl;
cin>>N;


// include the ghost inside initilization

double a[4]={0.0,0.0,0.0,0.0};

uint size;

double *d=nullptr;
double *X=nullptr;
double *Y=nullptr;

#if(!LID)
double Xa=-.5,Xb=1.0;
double Ya=-.5,Yb=0.5;
#else
double Xa=0.0,Xb=1.0;
double Ya=0.0,Yb=1.0;
#endif


const int level=2;



double val=0.0;
double x,y,z;
double dx=(Xb-Xa)/N;
double dy=(Yb-Ya)/N;
double dh=dx;
double C1=1./6.;
static const double dh2=dh*dh;

Q q(N,dx,dy);

q.initialize(N,a);

//q.wGetFace(1, &size, &d );

delete[] d;


double Re=20.0;
double uold;



int count=0;

// calculate average distance
double D[3];
count=0;

double fv=.1;

double res=1.0;

#if(!LID)
q.setExactBC(Xa,Xb,Ya,Yb);

//q.showExact(Xa,Ya);

//q.update();
//q.Grad();

q.start();

//q.getResTotal(Re);
//q.getResNorm(&res);
//q.debug(Xa,Ya);
//cout<<"Res "<<res<<endl;
 count=0;
//for(int i=0;i<1000000;i++)
while(res>1.e-12)
{
q.getRes(Re);
q.predict();
q.project();
q.setNeumanPressure();
//q.setExactBC(Xa,Xb,Ya,Yb);
q.correct();
q.update();

q.getResTotal(Re);
q.getResNorm(&res);
count++;
if(count % 60 ==0)
{
cout<<res<<endl;
}
}

#else
// solve lid driven Cavity


//q.showExact(Xa,Ya);
Re=1000.;

q.setBoundaryLidDrivenCavity();
//q.update();
//q.Grad();

q.start();

//q.getResTotal(Re);
//q.getResNorm(&res);
//q.debug(Xa,Ya);
//cout<<"Res "<<res<<endl;
 count=0;
//for(int i=0;i<1000000;i++)
#if(1)
double res_old=2.0;
while(res>1.e-12)
//for(int i=0;i<50000;i++)
{
q.getRes(Re);
q.predict();
q.setNeumanPressureLDC();
q.project();
q.setBoundaryLidDrivenCavity();
q.correct();
q.update();

//q.setBoundaryLidDrivenCavity();
res_old=res;

q.getResTotal(Re);
q.getResNorm(&res);
count++;
if(count % 100 ==0)
{
//cout<<count<<" "<<res<<endl;
cout<<res<<endl;
}
if(fabs(res_old-res)<1e-14)
{
//break;
}
res_old=res;
}
#endif


#endif












cout<<res<<endl;

/*
q.uSetFace(&fv);
//


//q.setSolid();

for(uint i=0;i<2000;i++)
//while(res>1e-12)
{
q.getRhsAB(Re);

q.poisson();

q.correct();

q.setNeumanPressure();

q.setSolid();

q.update(&res);

//cout<<i<<" "<<res<<endl;
cout<<" "<<res<<endl;

}
*/

#if(0)
Struct_2D(Xa,Xb,Ya,Yb, N,N,&X, &Y);
q.VTK_out(X,Y,N);

#else
q.Struct_2D_Ghost(Xa,Xb,Ya,Yb,&X,&Y);
q.VTK_out_with_ghost(X,Y);
#endif
};


void Struct_2D(double Xa,double Xb, double Ya, double Yb,int N, int M,double **X, double **Y)
{
	int i;
	
	double hx=N;
	double hy=N;
        
	double hz=N;
	double Xh=(Xb-Xa)/(hx);
	double Yh=(Yb-Ya)/(hy);

        cout<<"delx "<<Xh<<endl;       
 
	double Zh=0.01*Xh;
	
	(*X)=(double*)new double[N+1];    
	(*Y)=(double*)new double[N+1];	

        	
	for(i=0;i<N+1;i++)
	{
	(*X)[i]=Xa+Xh*i;		
	}

	for(i=0;i<N+1;i++)
	{
	(*Y)[i]=Ya+Yh*i;	
	}

}



double rhs(double &x, double &y, double &z)
{
static const double pi=atan(1.0)*4.0;
static const double pi_2=pi*pi;
return(-(3.0*pi_2*sin((pi*x)/2.0)*sin((pi*y)/2.0)*sin((pi*z)/2.0))/4.0);
}




