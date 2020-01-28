Coordinate_system *coord; 
Fields *fields;
Gas *gas;
double dt;
double c = 137.035999084;

TDSE::TDSE
// Constructor 
TDSE(Coordinate_system* coord, Fields* fields, Gas* gas, double dt);
// Destructor
~TDSE();

void timestep();

void dipole();
