#include "Coordinate_system.h"
#include <math.h>       /* floor */
#include <omp.h>

Coordinate_system::Coordinate_system(double z_max, double dz, double x_max, double dx)
{
  this->z_max = z_max;
  this->dz = dz;
  this->x_max = x_max;
  this->dx = dx;

  nz = 2*floor(z_max/dz)+1;
  z = new double[nz];
  for(int i = 0; i<nz;i++)
  {
    z[i] = i*dz-z_max;
  }

  nx = 2*floor(x_max/dx)+1;
  x = new double[nx];
  for(int i = 0; i<nx;i++)
  {
    x[i] = i*dx-x_max;
  }

}

Coordinate_system::~Coordinate_system()
{
  delete[] z;
  delete[] x;
}
