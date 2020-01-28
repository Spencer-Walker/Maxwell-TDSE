#include "Complex.h"

Complex::Complex()
{
}


Complex::Complex(double* z, int n)
{
  this->n = n;
  this->z = z;
}

Complex::~Complex()
{

}

void Complex::operator += (Complex obj)
{
  for(int i = 0;i<2*n;i++)
  {
    z[i]+= obj.z[i];
  }
}

void Complex::operator += (double* x)
{
  for(int i = 0;i<n;i++)
  {
    z[2*i] += x[i];
  }
}

void Complex::operator %= (double* y)
{
  for(int i = 0;i<n;i++)
  {
    z[2*i+1] += y[i];
  }
}

void Complex::operator -= (Complex obj)
{
  for(int i = 0;i<2*n;i++)
  {
    z[i]-= obj.z[i];
  }
}

void Complex::operator *= (Complex obj)
{
  if (obj.n == 1)
  {
    for(int i = 0;i<n;i++)
    {
      double tmp = z[2*i];
      z[2*i] = tmp*obj.z[0] - z[2*i+1]*obj.z[1];
      z[2*i+1] = tmp*obj.z[1] + z[2*i+1]*obj.z[0];
    }
  }
  else if(obj.n == n)
  { 
    for(int i = 0;i<n;i++)
    {
      double tmp = z[2*i];
      z[2*i] = tmp*obj.z[2*i] - z[2*i+1]*obj.z[2*i+1];
      z[2*i+1] = tmp*obj.z[2*i+1] + z[2*i+1]*obj.z[2*i];
    }
  }
}

void Complex::operator *= (double* x)
{
  for(int i = 0;i<n;i++)
  {
    z[2*i] *= x[i];
    z[2*i+1] *= x[i];
  }
}

void Complex::operator *= (double x)
{
  for(int i = 0;i<n;i++)
  {
    z[2*i] *= x;
    z[2*i+1] *= x;
  }
}

void Complex::operator /= (Complex obj)
{
   if (obj.n == 1)
  {
    double tmp2 = obj.z[0]*obj.z[0]+obj.z[1]*obj.z[1];
    for(int i = 0;i<n;i++)
    {
      double tmp1 = z[2*i];
      z[2*i] = (tmp1*obj.z[0] + z[2*i+1]*obj.z[1])/tmp2;
      z[2*i+1] = (-tmp1*obj.z[1] + z[2*i+1]*obj.z[0])/tmp2;
    }
  }
  else if(obj.n == n)
  { 
    for(int i = 0;i<n;i++)
    {
      double tmp1 = z[2*i];
      double tmp2 = obj.z[2*i]*obj.z[2*i]+obj.z[2*i+1]*obj.z[2*i+1];
      z[2*i] = (tmp1*obj.z[2*i] + z[2*i+1]*obj.z[2*i+1])/tmp2;
      z[2*i+1] = (-tmp1*obj.z[2*i+1] + z[2*i+1]*obj.z[2*i])/tmp2;
    }
  }
}

void Complex::operator /= (double* x)
{
  for(int i = 0;i<n;i++)
  {
    z[2*i] /= x[i];
    z[2*i+1] /= x[i];
  }
}

void Complex::operator /= (double x)
{
  for(int i = 0;i<n;i++)
  {
    z[2*i] /= x;
    z[2*i+1] /= x;
  }
}

double Complex::zr(int i)
{
  return z[2*i];
}

void Complex::zr(double x,int i)
{
  z[2*i] = x;
}

double Complex::zi(int i)
{
  return z[2*i+1];
}

void Complex::zi(double y, int i)
{
  z[2*i+1] = y;
}

double Complex::r2(int i)
{
  return z[2*i]*z[2*i] + z[2*i+1]*z[2*i+1];
}

void Complex::set_z(double* z, int n)
{
  this->z = z;
  this->n = n;
}