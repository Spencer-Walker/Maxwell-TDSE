#include "Fields.h"
#include <omp.h>
Fields::Fields(Coordinate_system* coord, double z0, double E0, double fwhm, double dt)
{
  this->coord = coord;
  double tn = 0.0;
  double tnm1 = -dt;
  double tnm2 = -2.0*dt;
  En   = new double[coord->nz];
  Enm1 = new double[coord->nz];
  Enm2 = new double[coord->nz];
  D2Pn = new double[coord->nz];
  Pn   = new double[coord->nz];

  for(int i = 0; i<coord->nz; i++)
  {
    En[i]   = gaussian_beam( coord->z[i]-z0, tn, E0, z0, fwhm);
    Enm1[i] = gaussian_beam( coord->z[i]-z0, tnm1, E0, z0, fwhm);
    Enm2[i] = gaussian_beam( coord->z[i]-z0, tnm2, E0, z0, fwhm);
    D2Pn[i] = 0.0;
    Pn[i] = 0.0;
  }

}

Fields::~Fields()
{
  delete[] En;
  delete[] Enm1;
  delete[] Enm2;
  delete[] D2Pn;
  delete[] Pn;
}

double Fields::gaussian_beam(double z, double t, double E0, double z0, double fwhm)
{
  return E0*pow(2.0,-4.0*pow(((z-c*t)/fwhm),2))*cos(z-c*t);
}
