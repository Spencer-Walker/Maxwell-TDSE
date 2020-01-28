#include "Gas.h"
#include <iostream>
#include <omp.h>
Gas::Gas(Coordinate_system* coord, double zs_max, double dzs, double n0, double fwhm)
{
  this->coord = coord;
  this->zs_max = zs_max;
  this->dzs = dzs;
  this->n0 = n0;
  this->fwhm = fwhm; 

  ns = 2*floor(zs_max/dzs)+1;
  zs = new double[ns];
  #pragma omp parallel for
  for(int i = 0; i<ns;i++)
  {
    zs[i] = i*dzs-zs_max;
  }

  density = new double[coord->nz];
  #pragma omp parallel for
  for(int i = 0; i<coord->nz; i++)
  {
    density[i] = n0*pow(2.0,-4.0*pow((coord->z[i]/fwhm),2));
  }


}

Gas::~Gas()
{
  delete[] zs;
  delete[] density;
}
