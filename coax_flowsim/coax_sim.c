#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"
#include "reduced.h"

//Define dimless parameters
#define Re_g 2666.021505
#define Re_l 49.588
#define We 3.031640962
#define Fr 214.5104818
#define D_rel 0.732484076
#define Rho_rel 0.001225
#define U_rel 4.356435644
#define domsize_coeff 10
//define repeating variables
#define U_g 1.
#define D_i 1.
#define Rho_g 1.

int maxlevel=10;					//max level of refinement
double uemax=0.01;					//max velocity error 
double femax=0.001; 				//max volume fraction error

/*
uemax and femax are important when choosing the computational domain size as
these two parameters determine when grid refinement is applied.
Generally, larger computational domain requires a smaller error threshold.
*/

/* -------------------imposing boundary conditions------------------------------*/
scalar f0[];
u.n[left] = dirichlet(U_g+((U_g/U_rel)-U_g)*((y-D_i)>0 ? 1. : 0.)-(U_g/U_rel)*((y-(D_i/D_rel))>0 ? 1. : 0.));
									//applying step velocity profile at inlet	
u.t[left] = dirichlet(0.); 									

#if dimension>2
u.r[left] = dirichlet(0.); 			//applying bc of zero radial velocity at left boundary
#endif


p[left] = neumann(0.); 				//applying bc of zero pressure gradient at left boundary
f[left] = f0[];

u.n[right] = neumann(0.); 			//applying bc of zero velocity gradient at right boundary
p[right] = dirichlet(0.);			//applying bc of zero pressure at right boundary (?) 

/* ----------------------main program---------------------------------*/
int main(int argc, char * argv[]){
	if (argc>1)
		maxlevel = atoi (argv[1]);	//optional cl arguments: maximum level
	if (argc >2)
		uemax = atof (argv[2]);		//optional cl arguments: error threshold

	init_grid(256); 				//discretise initial domain with 256^2 grid points
	size(domsize_coeff*D_i);
	DT = 0.1;

	rho1=Rho_g/Rho_rel,rho2=Rho_g;	//Note  subscript 2 is gas, 1 is water
	mu1=D_i*rho2*U_g/Re_l;
	mu2=D_i*rho2*U_g/Re_g;
	f.sigma=D_i*sq(U_g)*rho2/We;

	G.x = sq(U_g)/(D_i*Fr);
	run();

}

/* -------------------setting initial conditions -------------------------------*/
event init (t=0){
	if(!restore (file="restart")){

		/*initial grid refinement in a cylinder at:
			99% of the length of computational domain 
			and 2 times the initial outer diameter*/
		refine (x<(domsize_coeff*0.99)*D_i && y<2.*D_i/D_rel && level <maxlevel);
		
		/*initialise annular liquid sheet volume fraction*/
		fraction (f0, difference((D_i/D_rel - y),(D_i - y)));
		f0.refine=f0.prolongation=fraction_refine;
		restriction ({f0});
		
		/*change the shape of initial jet to a half-circular shape
		this is done to reduce effects of surface tension on the initial jet, 
		then apply the initial jet velocity*/
		foreach(){
			f[]=f0[]*(0<x && D_i<y && y<D_i/D_rel && x<sqrt(sq((D_i/D_rel-D_i)/2.)-sq(y-D_i-(D_i/D_rel-D_i)/2.)));
			u.x[]=f[]*U_g/U_rel;
		}
	}
}

/* -------------------Output -------------------*/
event logfile (i++){
	if (i==0)
			fprintf (stderr, "t dt mgp.i mgpf.p mgu.i grid->tn perf.t perf.speed\n");
	fprintf (stderr, "%g %g %d %d %d %ld %g %g\n",t,dt,mgp.i,mgpf.i,mgu.i,grid->tn,perf.t,perf.speed);
	
}

event movie (t += 1e-2){
	view (tx = -0.5);
	clear();
	draw_vof("f");
	//squares ("u.x", linear = true);
	box();
	mirror(n={0,1}){
		draw_vof("f");
		//squares ("u.x", linear = true);
		box();	
	}
	save("movie.mp4");
}

/*Labels and outputs volume fraction at specified timestep as binary files
See accompanying python scripts to read binary file*/
event field_binout (t+=0.1){
	char filename[256];
	sprintf(filename,"frac-%f.bin",t);
	FILE * fp_f=fopen(filename,"w");
	output_matrix(f, fp_f, linear=true);
	fclose(fp_f);
}


event end (t=0.;t+=0.1;t<=15.){
}

event adapt(i++){
	adapt_wavelet({f,u}, (double[]){femax,uemax,uemax,uemax},maxlevel)

	/*THis is a stop gap solution to prevent backflow during simulation
	by coarsening the mesh to the coarsest setting after 99% of computational domain length.
	However, backflow are usually caused by incorrect femax and uemax values or other parameters,
	and volume might not be conserved 
	so this measure should only be used when necessary*/

	//unrefine (x>(domsize_coeff*0.99)*D_i);	

}