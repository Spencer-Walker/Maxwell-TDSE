// mpirun -n 1 -genv OMP_NUM_THREADS=4 -genv I_MPI_PIN_DOMAIN=omp ./Maxwell-TDSE
#include <omp.h>
#include "mpi.h"
#include <math.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <fstream>
#include "Maxwell.h"
#include "Fields.h"
#include "Gas.h"
#include "Linear_algebra.h"
#include "Complex.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include "TDSE.h"
using namespace Eigen;

int main(int argc, char *argv[]) {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int iam = 0, np = 1;

	// Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

	// Test parallel configuration 
  #pragma omp parallel default(shared) private(iam)
  {
    np = omp_get_num_threads();
    iam = omp_get_thread_num();
    printf("Hello from thread %d out of %d from process %d out of %d on %s\n",
            iam, np, rank, numprocs, processor_name);
  }

	// Maxwell Geometry (z) and TDSE Geometry (x) (Everythin 1D for now)	
  double z_max = 150;        // Maxwell box size [-z_max,z_max]
  double dz = 0.5;           // Maxwell grid step
  double r = 1.0;			       // Maxwell CFL
  double dt = dz*r/137.0;    // Time step (both share for now)
  int nt = floor(1.0/dt);	 // Number of time steps
  double x_max = 50;				 // TDSE box size [-x_max,x_max]
  double dx = 0.1;					 // TDSE grid step 
  Coordinate_system coord(z_max, dz,x_max, dx); // Initialization

	// Output TDSE grid.
	std::ofstream x_file ("x.out");
	for(int ix = 0; ix < coord.nx; ix++)
	{
		x_file << coord.x[ix] << " , ";
	}
	x_file.close();

	// Output Maxwell grid.
	std::ofstream z_file ("z.out");
	for(int iz = 0; iz < coord.nz; iz++)
	{
		z_file << coord.z[iz] << " , ";
	}
	z_file.close();

	// Initialize the gausian pulse
  double z0 = -75.0;                 // Initial pulse location
  double E0 = 0.05;                 // Initial Electric field strength
  double fwhm_field = 20.0;        // Full width half max of the field.
  Fields fields(&coord, z0, E0, fwhm_field, dt); // Initialization

	// Output initial electric field.
	std::ofstream E0_file ("E0.out");
	for(int iz = 0; iz < coord.nz; iz++)
  {
    E0_file << fields.En[iz] << " , ";
  }
	E0_file.close();

	// Distribute the TDSE points.
  double zs_max = 40;      // Sample box size [-zs_max,zs_max]
  double dzs = 0.1;			   // Sample grid step
  double n0 = 50.0;			   // Gas density peak
  double fwhm_gas = 20.0; // Full width half max gas density
  Gas gas(&coord, zs_max, dzs, n0, fwhm_gas); // Initialization

	// Output sample points zs
	std::ofstream zs_file ("zs.out");
	for(int is = 0; is < gas.ns; is++)
	{
		zs_file << gas.zs[is] << " , ";
	}	
	zs_file.close();

	// Output density (at all z)
	std::ofstream n_file ("n.out");
	for(int iz = 0; iz < coord.nz; iz++)
  {
   	n_file << gas.density[iz] << " , ";
  }
	n_file.close();

  // Initialize linear algebra package 
  Linear_algebra lin = Linear_algebra();
	
	// Initialize maxwell solver
	Maxwell maxwell(&coord,&fields, &gas, dt);
  
	// Initialize the TDSE solver 
	TDSE tdse = TDSE(&coord, &fields, &gas, dt);
     
	// Allocate wavefunctions for the samples
  double** psi = new double*[gas.ns];  
  Complex PSI[gas.ns];
  for(int is = 0; is<gas.ns; is++)
  {
    psi[is] = new double[2*coord.nx]; // 2x since complex
    PSI[is].set_z(psi[is],coord.nx); 
  }

	// Each thread gets its own "potential"
  double** V = new double*[np];
  Complex VV[np];
	
	// B is needed for each linear solve.
  double** b = new double*[np];
	Complex B[np];

	// RHS is the old wavefunction and then the rhs of CN. (more workspace)
	double** rhs = new double*[np];
  Complex RHS[np];
	
	// Allocate everything that each rank needs.
	#pragma omp parallel default(shared) private(iam)
  {
		iam = omp_get_thread_num();
    V[iam] = new double[2*coord.nx];
    VV[iam].set_z(V[iam],coord.nx);
    b[iam] = new double[2*coord.nx];
    B[iam].set_z(b[iam],coord.nx);
  	rhs[iam] = new double[2*coord.nx];
    RHS[iam].set_z(rhs[iam],coord.nx);
	}

	//////////////// PUT THIS IN A CLASS ////////////////////////////////////////////////
  VectorXd s(coord.nx-1);
  for(int ix=0; ix<coord.nx-1; ix++)
  {
    s[ix] = -0.5/pow(coord.dx,2);
  }

  VectorXd d(coord.nx);
  for(int ix=0; ix<coord.nx; ix++)
  {
    d[ix] = 1/pow(coord.dx,2) + -1/sqrt(2+coord.x[ix]*coord.x[ix]);
  }

  SelfAdjointEigenSolver<MatrixXd>* eig = new SelfAdjointEigenSolver<MatrixXd>;
  
  eig->computeFromTridiagonal(d,s);

  VectorXd E = eig->eigenvalues();
  MatrixXcd states = eig->eigenvectors();

	std::ofstream states_file ("states.out");
	for(int state = 0; state<coord.nx; state++)
	{
		states_file << states.col(state).real().transpose() << "\n";
	}
	states_file.close();

  double* ground_state = new double[2*coord.nx];
  for(int ix = 0; ix<coord.nx; ix++)
  {
    ground_state[2*ix] = states.col(0)[ix].real();
    ground_state[2*ix+1] = 0;
  }

	
	for(int is =0; is<gas.ns; is++)
	{
		iam = omp_get_thread_num();
		for(int ix =0; ix<2*coord.nx; ix++)
		{
			PSI[is].z[ix] = ground_state[ix];
		}
	}

  states.resize(0,0);
  E.resize(0);
  delete eig;
  d.resize(0);
  s.resize(0);
//////////////////////////////////////////////////////////////////////////////////////

	// Electric feild file (all timesteps)
	std::ofstream E_file ("E.out");
	int iz = 0;	
	// Integrate the Maxwell-TDSE in time (put in class?)
  for(int it=0;it<nt;it++)
  { 
		// Print the time-step
    std::cout << "t = " << it*dt << std::endl; 
    // Integrate Maxwell one time-step
		maxwell.timestep();

		// Solve the TDSE for each sample in parallel
    #pragma omp parallel for private(iz,iam,np)
    for(int is = 0; is<gas.ns; is++)
    {
			// Which thread? 
      iam = omp_get_thread_num();
      // Map the sample onto a gridpoint.
			iz = floor((is*gas.dzs + coord.z_max - gas.zs_max)/coord.dz);    
			// Integrate each TDSE in time.
			tdse.timestep(&VV[iam], &B[iam], &PSI[is], &RHS[iam], &lin, dt, fields.En[iz]);
      // Set the Polarization density acceleration. 
			fields.D2Pn[iz] = gas.density[iz]*fields.En[iz];	
		}
		for(iz = 0; iz<coord.nz; iz++)
		{
			E_file << fields.En[iz] << " , ";
		}
		E_file << "\n";
	}
	E_file.close();

  std::cout << "DONE" << std::endl;
  std::cout << "nz = " << coord.nz << std::endl;
  std::cout << "nx = " << coord.nx << std::endl;
  std::cout << "ns = " << gas.ns << std::endl;
  std::cout << "nt = " << nt << std::endl;
  
  MPI_Finalize();
}
