#ifndef Coordinate_system_H
#define Coordinate_system_H
class Coordinate_system
{
  public:
    double* z;
    double z_max;
    double dz;
    int nz;
    double* x;
    double x_max;
    double dx;
    int nx;

    // Constructor
    Coordinate_system(double z_max, double dz, double x_max, double dx);

    // Destructor
    ~Coordinate_system();
};
#endif