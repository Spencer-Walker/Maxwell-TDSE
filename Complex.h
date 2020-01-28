#include <cmath>
#include <fstream>
#include <omp.h>
#include <iostream>

#ifndef Complex_H
#define Complex_H
class Complex
{
  public:
    double* z;
    int n;
    // Constructor
    Complex(); 
    Complex(double* z, int n);
    // Destructor
    ~Complex();
    void operator += (Complex obj);
    void operator += (double* x);
    void operator %= (double* y);
    void operator -= (Complex obj);
    void operator *= (Complex obj);
    void operator *= (double* x);
    void operator *= (double x);
    void operator /= (Complex obj);
    void operator /= (double* x);
    void operator /= (double x);
    double zr(int i);
    void zr(double x, int i);
    double zi(int i);
    void zi(double y,int i);
    double r2(int i);
    void set_z(double* z, int n);
};
#endif

