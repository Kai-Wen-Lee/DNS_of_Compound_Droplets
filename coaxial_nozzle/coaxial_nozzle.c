#inlcude "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

#define radius_outer 12.56/1000.
#define radius_inner 9.02/1000.

#define Re 5800			//to be changed to real value
#define SIGMA 3e-5		//to be changed to real value

int maxlevel = 10;
double uemax = 0.1;

//Define aux vol fraction field
scalar vf_0[];

//Set boundary conditions
u.n[left] = dirichlet(vf_0[]*1.);		//velocity tbc
u.t[left] = dirichlet(0);

p[left] = neumann(0);
f[left] = vf_0[];						//what does this line mean physically

u.n[right] = neumann(0);
p[right] = dirichlet(0);

int main (){
	init_grid (128);
	size (0.5);
	
	//setting density of each phase
	rho_l = 1000.;			//water
	rho_g = 1.225;			//air
	
	//setting viscosity
	mu_l = 2.*radius/Re*rho_l;
	mu_g = 2.*radius/Re*rho_g;
	
	//setting surface tension coefficient'
	f.sigma = SIGMA;
	
	//run
	run();
}

//set initial conditions
