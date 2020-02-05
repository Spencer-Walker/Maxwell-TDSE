#include "Linear_algebra.h"


// Constructor 
Linear_algebra::Linear_algebra()
{

}

// Destructor
Linear_algebra::~Linear_algebra()
{

}

void Linear_algebra::real_thomas_algorithm(double *a, double *b, double *c, double *d, double *x, int n)
{
	double w;
	for(int i=1;i<n;i++)
	{
		w = a[i-1]/b[i-1];
		b[i] = b[i]-w*c[i-1];
		d[i] = d[i]-w*d[i-1];
	}
	x[n-1] = d[n-1]/b[n-1];
	for(int i=n-2;i>=0;i--)
	{
		x[i] = (d[i]-c[i]*x[i+1])/b[i];
	}
}

void Linear_algebra::complex_thomas_algorithm(Complex* A, Complex* B, Complex* C, Complex* D, Complex* X, int n)
{
	double* a = A->z;
	double* b = B->z;
	double* c = C->z;
	double* d = D->z;
	double* x = X->z;
	double wr;
	double wi;

	double tmp;

	for(int i=1;i<n;i++)
	{
		tmp = B->r2(i-1);
		wr = (a[2*(i-1)]*b[2*(i-1)] + a[2*(i-1)+1]*b[2*(i-1)+1])/tmp;
		wi = (a[2*(i-1)+1]*b[2*(i-1)] - a[2*(i-1)]*b[2*(i-1)+1])/tmp;

		b[2*i]   = b[2*i]-(wr*c[2*(i-1)] - wi*c[2*(i-1)+1]);
		b[2*i+1] = b[2*i+1]-(wr*c[2*(i-1)+1] + wi*c[2*(i-1)]);

		d[2*i]   = d[2*i]  -(wr*d[2*(i-1)] - wi*d[2*(i-1)+1]);
		d[2*i+1] = d[2*i+1]-(wr*d[2*(i-1)+1] + wi*d[2*(i-1)]);
	}
	tmp = B->r2(n-1);
	x[2*(n-1)]   = (d[2*(n-1)]*b[2*(n-1)] + d[2*(n-1)+1]*b[2*(n-1)+1])/tmp;
	x[2*(n-1)+1] = (d[2*(n-1)+1]*b[2*(n-1)] - d[2*(n-1)]*b[2*(n-1)+1])/tmp;
	double tmpr;
	double tmpi;
	for(int i=n-2;i>=0;i--)
	{
		tmpr = d[2*i]-(c[2*i]*x[2*(i+1)]-c[2*i+1]*x[2*(i+1)+1]);
		tmpi = d[2*i+1]-(c[2*i]*x[2*(i+1)+1]+c[2*i+1]*x[2*(i+1)]);
		tmp  = B->r2(i);
		x[2*i]   = (tmpr*b[2*i] + tmpi*b[2*i+1])/tmp;
		x[2*i+1] = (tmpi*b[2*i] - tmpr*b[2*i+1])/tmp;
	}
		
}

void Linear_algebra::mat_vec_real_tri_mul(double* a, double* b, double* c, double* x, int n)
{
	double t0;
	double t1;

	t0 = x[0];
	x[0] = b[0]*t0 + c[0]*x[1];
	for(int i = 1; i<n-1; i++)
	{
		t1 = x[i];
		x[i] = a[i-1]*t0 + b[i]*t1 + c[i]*x[i+1];
		t0 = t1;
	}
	t1 = x[n-1];
	x[n-1] = a[n-2]*t0 + b[n-1]*t1;

}

void Linear_algebra::mat_vec_complex_tri_mul(Complex* A, Complex* B, Complex* C, Complex* X, int n)
{
	double t0r;
	double t0i;
	double t1r;
	double t1i;


	double* a = A->z;
	double* b = B->z;
	double* c = C->z;
	double* x = X->z;

	t0r = x[0];
	t0i = x[1];

	x[0] = b[0]*t0r - b[1]*t0i + c[0]*x[2] - c[1]*x[3];
	x[1] = b[0]*t0i + b[1]*t0r + c[0]*x[3] + c[1]*x[2]; 

	for(int i = 1; i<n-1; i++)
	{
		t1r = x[2*i];
		t1i = x[2*i+1];
		x[2*i] = a[2*(i-1)]*t0r - a[2*(i-1)+1]*t0i + b[2*i]*t1r - b[2*i+1]*t1i + c[2*i]*x[2*(i+1)] - c[2*i+1]*x[2*(i+1)+1];
		x[2*i+1] = a[2*(i-1)]*t0i + a[2*(i-1)+1]*t0r + b[2*i]*t1i + b[2*i+1]*t1r + c[2*i]*x[2*(i+1)+1] + c[2*i+1]*x[2*(i+1)];

		t0r = t1r;
		t0i = t1i;
	}

	t1r = x[2*(n-1)];
	t1i = x[2*(n-1)+1];
	x[2*(n-1)]   = a[2*(n-2)]*t0r - a[2*(n-2)+1]*t0i + b[2*(n-1)]*t1r - b[2*(n-1)+1]*t1i;
	x[2*(n-1)+1] = a[2*(n-2)]*t0i + a[2*(n-2)+1]*t0r + b[2*(n-1)]*t1i + b[2*(n-1)+1]*t1r;

}

void Linear_algebra::cn_rhs(double dt,double dx, Complex* VV, Complex* PSI, int n)
{
	double t0r;
	double t0i;
	double t1r;
	double t1i;

	double a =  (dt/(4*pow(dx,2)));
	double c = 1/pow(dx,2);
	double* V = VV->z;
	double* psi = PSI->z;

	t0r = psi[0];
	t0i = psi[1];

	psi[0] = (1+dt*V[1]/2)*t0r - (-dt*(c + V[0])/2 )*t0i - a*psi[3];
	psi[1] = (1+dt*V[1]/2)*t0i + (-dt*(c + V[0])/2 )*t0r + a*psi[2]; 
	for(int i = 1; i<n-1; i++)
	{
		t1r = psi[2*i];
		t1i = psi[2*i+1];
		psi[2*i] =   - a*t0i + (1+dt*V[2*i+1]/2)*t1r - (-dt*(c + V[2*(n-1)])/2 )*t1i - a*psi[2*(i+1)+1];
		psi[2*i+1] = + a*t0r + (1+dt*V[2*i+1]/2)*t1i + (-dt*(c + V[2*(n-1)])/2 )*t1r + a*psi[2*(i+1)];
		t0r = t1r;
		t0i = t1i;
	}

	t1r = psi[2*(n-1)];
	t1i = psi[2*(n-1)+1];
	psi[2*(n-1)]   = - a*t0i + (1+dt*V[2*(n-1)+1]/2)*t1r - (-dt*(c + V[2*(n-1)])/2 )*t1i;
	psi[2*(n-1)+1] = + a*t0r + (1+dt*V[2*(n-1)+1]/2)*t1i + (-dt*(c + V[2*(n-1)])/2 )*t1r;

}

void Linear_algebra::cn_lhs(double dt,double dx, Complex* VV, Complex* B, Complex* PSI, Complex* RHS, int n)
{
	double wr;
	double wi;
	double* V = VV->z;
	double* b = B->z;
	double* psi = PSI->z;
	double* rhs = RHS->z;

	double tmp;
	double a =  -(dt/(4*pow(dx,2)));
	double c = 1/pow(dx,2);
	double dt2 = dt/2;
	for(int i=0;i<n;i++)
	{
		b[2*i] = (1-dt2*V[2*i+1]);
		b[2*i+1] = dt2*(c + V[2*i]);
	}

	for(int i=1;i<n;i++)
	{
		tmp = B->r2(i-1);
		wr = (a*b[2*(i-1)+1])/tmp;
		wi = (a*b[2*(i-1)])/tmp;

		b[2*i]   = b[2*i]+wi*a;
		b[2*i+1] = b[2*i+1]-wr*a;

		rhs[2*i]   = rhs[2*i]-(wr*rhs[2*(i-1)] - wi*rhs[2*(i-1)+1]);
		rhs[2*i+1] = rhs[2*i+1]-(wr*rhs[2*(i-1)+1] + wi*rhs[2*(i-1)]);
	}
	tmp = B->r2(n-1);
	psi[2*(n-1)]   = (rhs[2*(n-1)]*b[2*(n-1)] + rhs[2*(n-1)+1]*b[2*(n-1)+1])/tmp;
	psi[2*(n-1)+1] = (rhs[2*(n-1)+1]*b[2*(n-1)] - rhs[2*(n-1)]*b[2*(n-1)+1])/tmp;
	double tmpr;
	double tmpi;
	for(int i=n-2;i>=0;i--)
	{
		tmpr = rhs[2*i] + a*psi[2*(i+1)+1];
		tmpi = rhs[2*i+1] - a*psi[2*(i+1)];
		tmp  = B->r2(i);
		psi[2*i]   = (tmpr*b[2*i] + tmpi*b[2*i+1])/tmp;
		psi[2*i+1] = (tmpi*b[2*i] - tmpr*b[2*i+1])/tmp;
	}
}
