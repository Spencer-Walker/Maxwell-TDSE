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

using namespace Eigen;

int main(int argc, char *argv[]) {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int iam = 0, np = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  #pragma omp parallel default(shared) private(iam)
  {
    np = omp_get_num_threads();
    iam = omp_get_thread_num();
    printf("Hello from thread %d out of %d from process %d out of %d on %s\n",
            iam, np, rank, numprocs, processor_name);
  }

	
  double z_max = 150;
  double dz = 0.1;
  double r = 1.0;
  double dt = dz*r/137.0;
  int nt = floor(0.001/dt);
  double x_max = 50;
  double dx = 0.1;
  Coordinate_system coord(z_max, dz,x_max, dx);

  double z0 = 0.0;
  double E0 = 1.0;
  double fwhm_field = 20.0;
  Fields fields(&coord, z0, E0, fwhm_field, dt);

  double zs_max = 50;
  double dzs = 0.1;
  double n0 = 1.0;
  double fwhm_gas = 100.0;
  Gas gas(&coord, zs_max, dzs, n0, fwhm_gas);

  double** psi = new double*[gas.ns];
  Complex PSI[gas.ns];
  for(int is = 0; is<gas.ns; is++)
  {
    psi[is] = new double[2*coord.nx];
    PSI[is].set_z(psi[is],coord.nx);
  }

  double** V = new double*[np];
  Complex VV[np];
  for(int ip = 0; ip<np; ip++)
  {
    V[ip] = new double[2*coord.nx];
    VV[ip].set_z(V[ip],coord.nx);
  }

  double** b = new double*[np];
  Complex B[np];
  for(int ip = 0; ip<np; ip++)
  {
    b[ip] = new double[2*coord.nx];
    B[ip].set_z(b[ip],coord.nx);
  }

  double** rhs = new double*[np];
  Complex RHS[np];
  for(int ip = 0; ip<np; ip++)
  {
    rhs[ip] = new double[2*coord.nx];
    RHS[ip].set_z(rhs[ip],coord.nx);
  }

  Linear_algebra lin = Linear_algebra();
  Maxwell maxwell(&coord,&fields, &gas, dt);
	
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

  double* ground_state = new double[2*coord.nx];
  for(int ix = 0; ix<coord.nx; ix++)
  {
    ground_state[2*ix] = states.col(0)[ix].real();
    ground_state[2*ix+1] = 0;
  }

  double norm = 0;
	int ix = 0;
	
	for(int is =0; is<gas.ns; is++)
	{
    norm = 0;
		iam = omp_get_thread_num();
		for(int ix =0; ix<2*coord.nx; ix++)
		{
			PSI[is].z[ix] = ground_state[ix];
      norm += PSI[iam].z[ix]*PSI[iam].z[ix];
		}
    printf("norm =  %4.2f \n", norm);
	}

  states.resize(0,0);
  E.resize(0);
  delete eig;
  d.resize(0);
  s.resize(0);
	double* tmp;
 	ix = 0;
	int iz = 0;
  for(int it=0;it<nt;it++)
  { 
    std::cout << "it = " << it << std::endl;
    maxwell.timestep();
    #pragma omp parallel for private(ix,iz,iam, np, tmp,norm)
    for(int is = 0; is<gas.ns; is++)
    {
      iam = omp_get_thread_num();
			norm = 0;
      for(ix = 0; ix<coord.nx; ix++)
			{ 
        iz = floor((is*gas.dzs + coord.z_max - gas.zs_max)/coord.dz);
        VV[iam].z[2*ix] = -1/sqrt(2+coord.x[ix]*coord.x[ix])+fields.En[iz]*coord.x[ix];
        VV[iam].z[2*ix+1] = 0;
				norm += PSI[iam].z[2*ix]*PSI[iam].z[2*ix]+PSI[iam].z[2*ix+1]*PSI[iam].z[2*ix+1];
      }
			std::cout << norm << std::endl;
      tmp = RHS[iam].z;
      RHS[iam].z = PSI[is].z;
      PSI[is].z = tmp;
      lin.cn_rhs(dt,dx,&VV[iam], &RHS[iam], coord.nx);
      lin.cn_lhs(dt,dx,&VV[iam], &B[iam], &PSI[is], &RHS[iam], coord.nx);
      iz = floor((is*gas.dzs + coord.z_max - gas.zs_max)/coord.dz);
      fields.D2Pn[iz] = gas.density[iz]*fields.En[iz];	
		}
	}
  std::ofstream myfile ("example.txt");
  for(int i = 0; i<coord.nz;i++)
  {
    myfile << fields.En[i] << std::endl;
  }

  std::cout << "DONE" << std::endl;
  std::cout << "nz = " << coord.nz << std::endl;
  std::cout << "nx = " << coord.nx << std::endl;
  std::cout << "ns = " << gas.ns << std::endl;
  std::cout << "nt = " << nt << std::endl;
  
  MPI_Finalize();
}
