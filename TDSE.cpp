#include "TDSE.h"

// Constructor 
TDSE::TDSE(Coordinate_system* coord, Fields* fields, Gas* gas, double dt)
{
	this->coord  = coord;
	this->fields = fields;
	this->gas = gas;
	this->dt = dt; 
}

TDSE::~TDSE()
{

}

void TDSE::timestep(Complex* VV, Complex* B, Complex* PSI, Complex* RHS, Linear_algebra* lin,
										double dt, double E)
{
	for(int ix = 0; ix<coord->nx; ix++)
  {
    VV->z[2*ix] = -1/sqrt(2+coord->x[ix]*coord->x[ix])+E*coord->x[ix];
    VV->z[2*ix+1] = 0;
  }
	tmp = RHS->z;
  RHS->z = PSI->z;
  PSI->z = tmp;
  lin->cn_rhs(dt,coord->dx, VV, RHS, coord->nx);
  lin->cn_lhs(dt,coord->dx, VV, B, PSI, RHS, coord->nx);

}

double TDSE::dipole(Complex* PSI)
{
	double d = 0;
	for(int ix = 0; ix < coord->nx; ix++)
	{
		d += PSI->r2(ix)*coord->x[ix];
	}
	return d;
}


