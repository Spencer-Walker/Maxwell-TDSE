#include "Coordinate_system.h"
#include "Fields.h"
#include "Gas.h"
#include <cmath>
#include <fstream>
#include <omp.h>

#ifndef TDSE_H
#define TDSE_H
class TDSE
{
  public:
    Coordinate_system *coord; 
    Fields *fields;
    Gas *gas;
    double dt;
    double c = 137.035999084;

    // Constructor 
    TDSE(Coordinate_system* coord, Fields* fields, Gas* gas, double dt);
    // Destructor
    ~TDSE();

    void timestep();

    double dipole();

    double dipole_velocity();

    double dipole_acceleration();n

};
#endif