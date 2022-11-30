#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

//Define dimless parameters
#define Re_g 323.5489
#define Re_l 8956.075
#define We 1.241522
#define Fr 17.81714
#define D_rel 0.718153
#define Rho_rel 0.001225
#define U_rel 7.06667


//define repeating variables
#define U_g 1.25
#define D_i 0.00902
#define Rho_g 1.225

int maxlevel=8;	//max level of refinement =10
double uemax=0.1;	//error threshold of velocity is 0.1

/* -------------------imposing boundary conditions------------------------------*/
scalar f0[];												 //set an aux volume fraction field (1 is inside cylinder, 0 is outside cylinder)
u.n[left] = dirichlet(f0[]);		
u.t[left] = dirichlet(0); 									

#if dimension>2
u.r[left] = dirichlet(0); 		//applying bc of zero radial velocity at left boundary
#endif


p[left] = neumann(0); 			//applying bc of zero pressure gradient at left boundary
f[left] = f0[];

u.n[right] = neumann(0); 			//applying bc of zero velocity gradient at right boundary
p[right] = dirichlet(0);			//applying bc of zero pressure at right boundary (?) 
/* ----------------------main program---------------------------------*/
int main(int argc, char * argv[]){
	if (argc>1)
		maxlevel = atoi (argv[1]);	//optional cl arguments: maximum level
	if (argc >2)
		uemax = atof (argv[2]);		//optional cl arguments: error threshold

	init_grid(64); //discretise initial domain with 64^3 grid points
	size(25*D_i);
	DT = 0.1;

	rho1=Rho_g/Rho_rel,rho2=Rho_g;	//Note  subscript 2 is gas, 1 is water
	mu1=(D_i/D_rel)*rho1*(U_g/U_rel)/Re_l;
	mu2=D_i*rho2*U_g/Re_g;
	f.sigma=D_i*sq(U_g);

	run();

}

/* -------------------setting initial conditions -------------------------------*/
event init (t=0){
	if(!restore (file="restart")){

		/*use a static refinement down to maxlevel in a cyl. 1.2 times longer than 
		the init jet and twice the radies*/
		refine (x<10*D_i && y<2.*D_i/D_rel && level <maxlevel);

		fraction (f0, difference((D_i/D_rel - y), (D_i - y)));
		f0.refine=f0.prolongation=fraction_refine;
		restriction ({f0});		//for BC on levels
		
		/* use this to define the init jet and its velocity
		foreach(){
			
			f[]=f0[]*10*D_i;
			u.x[]=f[]*U_g/U_rel;
		}*/
	}
}
/*
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
  av.x[] += f[]*Fr*sq(U_g)/D_i;
}*/

/* -------------------Output -------------------*/

event logfile (i++){
	int j=0;
	foreach_boundary(left){
		j++;
		fprintf(stderr,"%d	%d	%f\n",i ,j ,u.x[0,1]);
	}
	
}

event end (t=0.1;t+=0.1;t<=3){
}

event adapt(i++){
	adapt_wavelet({f,u}, (double[]){0.01,uemax,uemax,uemax},maxlevel);
}
