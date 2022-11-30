#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

//Define radius of jet, init jet len, Re and surface tention coeff

#define radius 1./12.
#define length 0.025
#define Re 5800
#define SIGMA 3e-5

int maxlevel=10;	//max level of refinement =10
double uemax=0.1;	//error threshold of velocity is 0.1

/* -------------------imposing boundary conditions------------------------------*/
scalar f0[];												 //set an aux volume fraction field (1 is inside cylinder, 0 is outside cylinder)
u.n[left] = dirichlet(f0[]*(1.+0.05*sin(10.*2.*pi*t)));		//applying oscilating inflow velocity on lhs
u.t[left] = dirichlet(0); 									//applying BC of zero tangential velocity at inflow
/*u.n[left] is the ghost value of the scalar field outside the left boundary
	here the dirichlet BC is applied and has a value of zero
*/

#if dimension>2
u.r[left] = dirichlet(0); 		//applying bc of zero radial velocity at left boundary
#endif


p[left] = neumann(0); 			//applying bc of zero pressure gradient at left boundary
f[left] = f0[];

u.n[right] = neumann(0); 			//applying bc of zero velocity gradient at right boundary
p[right] = dirichlet(0);			//applying bc of zero pressure at right boundary (?) 
/* ----------------------main program---------------------------------*/
int main(){
	run();

}

/* -------------------Output -------------------*/

event logfile (i++){
	foreach()
	fprintf("i:	", i, "\n", f0[]);
}