#include "Coordinate_system.h"
#include "Fields.h"
#include "Gas.h"
#include <cmath>
#include <fstream>
#include "Complex.h"
#include "Linear_algebra.h"
#include <omp.h>

#ifndef TDSE_H
#define TDSE_H
class TDSE
{
  public:
    Coordinate_system *coord; 
    Fields *fields;
    Gas *gas;
    double* tmp;
		double dt;
    double c = 137.035999084;

    // Constructor 
    TDSE(Coordinate_system* coord, Fields* fields, Gas* gas, double dt);
    // Destructor
    ~TDSE();

    void timestep(Complex* VV, Complex* B, Complex* PSI, Complex* RHS, Linear_algebra* lin,
                    double dt, double E);

		double dipole(Complex* PSI);

    //double dipole_velocity();

    //double dipole_acceleration();

};
#endif
