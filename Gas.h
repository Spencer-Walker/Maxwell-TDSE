#include "Coordinate_system.h"

#ifndef Gas_H
#define Gas_H

#include <math.h>
class Gas
{
  public:
    Coordinate_system* coord;
    double zs_max;
    int ns;
    double* zs;
    double dzs;
    double* density;
    double n0;
    double fwhm;

    // Constructor
    Gas(Coordinate_system* coord, double zs_max, double dzs, double n0, double fwhm);

    // Destructor
    ~Gas();

};
#endif