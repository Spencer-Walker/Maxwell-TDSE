#include "Coordinate_system.h"
#include "Fields.h"
#include "Gas.h"
#include <cmath>
#include <fstream>
#include <omp.h>

#ifndef Maxwell_H
#define Maxwell_H
class Maxwell
{
  public:
    Coordinate_system *coord; 
    Fields *fields;
    Gas *gas;
    double* Etemp;
    double dt;
    double c = 137.035999084;
    double r;
    // Constructor 
    Maxwell(Coordinate_system* coord, Fields* fields, Gas* gas, double dt);
    // Destructor
    ~Maxwell();

    void timestep();



};
#endif