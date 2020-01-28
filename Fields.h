#include "Coordinate_system.h"
#include <math.h>

#ifndef Fields_H
#define Fields_H
class Fields
{
  public:
    Coordinate_system* coord;
    double c = 137.035999084;
    double* Enm2;
    double* Enm1;
    double* En;
    double* D2Pn;
    double* Pn;

    // Constructor
    Fields(Coordinate_system* coord, double z0, double E0, double fwhm, double dt);
    // Destructor
    ~Fields();

    double gaussian_beam(double z, double t, double E0, double z0, double fwhm);

};
#endif