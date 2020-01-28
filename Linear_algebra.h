#include <cmath>
#include <fstream>
#include <omp.h>
#include <iostream>
#include "Complex.h"

#ifndef Linear_algebra_H
#define Linear_algebra_H
class Linear_algebra
{
  public:
    double *a;
    double *b;
    double *c;
    double *d;
    double *x;

    // Constructor 
    Linear_algebra();
    // Destructor
    ~Linear_algebra();

    void real_thomas_algorithm(double *a, double *b, double *c, double *d, double *x, int n);
    void complex_thomas_algorithm(Complex* A, Complex* B, Complex* C, Complex* D, Complex* X, int n);
    void mat_vec_real_tri_mul(double* a, double* b, double* c, double* x, int n);
    void mat_vec_complex_tri_mul(Complex* A, Complex* B, Complex* C, Complex* X, int n);
    void cn_rhs(double dt,double dx, Complex* VV, Complex* PSI, int n);
    void cn_lhs(double dt,double dx, Complex* VV, Complex* B, Complex* PSI, Complex* RHS, int n);

};
#endif