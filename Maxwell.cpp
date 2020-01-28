#include "Maxwell.h"
#include <omp.h>
Maxwell::Maxwell(Coordinate_system* coord, Fields* fields, Gas* gas, double dt)
{
  this->coord = coord;
  this->fields = fields;
  this->gas = gas; 
  this->dt = dt;
  r = c*dt/coord->dz;
}

Maxwell::~Maxwell()
{

}

void Maxwell::timestep()
{
  Etemp = fields->Enm2;
  fields->Enm2 = fields->Enm1;
  fields->Enm1 = fields->En;
  fields->En = Etemp; 

  fields->En[0] = pow(r,2.0)*(fields->Enm1[1]-2.0*fields->Enm1[0])
                 +2*fields->Enm1[0]-fields->Enm2[0];
  int iam = -1;
  #pragma omp parallel for
  for(int iz = 1; iz<coord->nz-1; iz++)
  {
    fields->En[iz] = pow(r,2.0)*(fields->Enm1[iz+1]-2.0*fields->Enm1[iz]
                    +fields->Enm1[iz-1])+2*fields->Enm1[iz]-fields->Enm2[iz] - 4*M_PI*fields->D2Pn[iz]/pow(c,2.0);
  }
  fields->En[coord->nz-1] = pow(r,2.0)*(-2.0*fields->Enm1[coord->nz-1]+fields->Enm1[coord->nz-2])
                    +2*fields->Enm1[coord->nz-1]-fields->Enm2[coord->nz-1];
}
