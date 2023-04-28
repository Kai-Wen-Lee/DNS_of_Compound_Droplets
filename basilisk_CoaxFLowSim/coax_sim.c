#include "axi.h"						//module for axisymmetric flow
#include "navier-stokes/centered.h"		//module for cell centered incompressible variable density navier stokes flow
#include "navier-stokes/perfs.h"		//module for logging
#include "two-phase.h"					//module for volume of fliuid method
#include "tension.h"					//module for surface tension calculation
#include "view.h"						//module for Basilisk view to output animations
#include "reduced.h"					//module for gravity for interfacial flow
#include "output_vtu_bin_foreach.h"		//custom module by cselcuk
#include <stdint.h>

//Provide dimensionless numbers here
#define Re_g 117.7367742
#define Re_l 2.189904
#define We 0.064541337
#define Fr 573.7619776
#define D_rel 0.5333333
#define Rho_rel 0.00129
#define U_rel 1.356249755
#define domsize_coeff 20
#define osc_amp -0.1
#define osc_freq 0.1
//Repeating variables set to unity
#define U_g 1.
#define D_i 1.
#define Rho_g 1.

//Initial jet profile eclipse shape factors
#define ell_sf_A 2.0
#define ell_sf_B 1.5

//Define simulation time
#define T_END 60.

int maxlevel=7;	//max level of refinement =10
double uemax=0.000001;	//error threshold of velocity is 0.01
double femax=0.0000001; //error threshold of volume fraction field is 0.001

//Boundary Conditions
scalar f0[];	//aux volume fraction field for liquid phase	
u.n[left] = dirichlet((y >= 0 && y <= D_i) ? U_g : (y > D_i && y <= D_i/D_rel) ? (U_g/U_rel) : 0);	//Fluid inflow velocity profile normal to boundary
u.t[left] = dirichlet((y >= 0 && y <= D_i) ? 0 : (y > D_i && y <= D_i/D_rel) ? osc_amp*U_g/U_rel*sin(osc_freq*2.*pi*t) : 0);	//Fluid inflow velocity tangential to boundary set to zero 										

#if dimension>2
u.r[left] = dirichlet(0.); 		//applying bc of zero radial velocity at left boundary
#endif

p[left] = neumann(0.);	//Zero pressure gradient at left boundary
f[left] = (y > D_i && y <= D_i/D_rel) ? f[] : 0;	//Define fluid volume fraction field at boundary

// Free flow condition for the rest of the boundaries
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
//Main program (simulation starts here)
int main(int argc, char * argv[]){
	if (argc>1)
		maxlevel = atoi (argv[1]);
	if (argc >2)
		uemax = atof (argv[2]);


	//Computational domain settings
	init_grid(512);				//Initial grid
	size(domsize_coeff*D_i);		//Computational domain size
	//DT = 0.1;						//enforcing a minimum timestep of 0.1

	//FLuid properties obtained from dimensionless numbers
	rho1=Rho_g/Rho_rel,rho2=Rho_g; 	//rho1 = gas density, rho2 = liquid density
	mu1=D_i*rho2*U_g/Re_l;			//mu1 = gas viscousity
	mu2=D_i*rho2*U_g/Re_g;			//mu2 = water viscousity
	f.sigma=D_i*sq(U_g)*rho2/We;	//f.sigma = liquid surface tension
	G.x = sq(U_g)/(D_i*Fr);			//G.x = gravitational acceleration

	//TOLERANCE = 1e-4;
	run();
}

//Initial conditions
event init (t=0){
	if(!restore (file="restart")){
		/*Initial refinement in the region defined as below up to the maximum refinement level*/
		//refine (x<(domsize_coeff*0.99)*D_i && y<2.*D_i/D_rel && level < maxlevel);
		refine (x < domsize_coeff*0.9 && sq(y) + sq(z) < 2.*sq(D_i/D_rel) && level < maxlevel);
		//fraction (f0, difference((D_i/D_rel - y),(D_i - y)));
		//fraction (f0, ((y >= 0 && y <= D_i) ? 0 : (y > D_i && y <= D_i/D_rel) ? x<=ell_sf_A*(sqrt(sq((D_i/D_rel-D_i)/2.)-sq(y-D_i-(D_i/D_rel-D_i)/2.))) : 0));
		//fraction (f0, difference((sq(D_i/D_rel) - sq(y) - sq(z)),(sq(D_i) - sq(y) - sq(z))));
		fraction(f0,sq(D_i/D_rel) - sq(y) - sq(z));
		//fraction (f0, ((y >= 0 && y <= D_i) ? x>=ell_sf_A*sqrt(sq(D_i)-sq(y)) : (y > D_i && y <= D_i/D_rel) ? x<=ell_sf_B*sqrt(sq(D_i/D_rel)-sq(y)) : 0));
		f0.refine=f0.prolongation=fraction_refine;
		restriction ({f0});
		
		foreach(){
			//f[]=f0[]*(0<=x && D_i<=y && y<=D_i/D_rel && x<=ell_sf_A*(sqrt(sq((D_i/D_rel-D_i)/2.)-sq(y-D_i-(D_i/D_rel-D_i)/2.)))); // semi circle on yaxis
			f[] = f0[]*(sq(D_i/D_rel) > sq(y) + sq(z) && sq(D_i) < sq(y) + sq(z) && x<=ell_sf_A*(sqrt(sq((D_i/D_rel-D_i)/2.)-sq(y-D_i-(D_i/D_rel-D_i)/2.))));
			//f[] = f0[]*(x < length); //simple rectangular jet length ((y > 0 && y < D_i) ? x>=sqrt(sq(D_i)-sq(y)) : (y > D_i && y < D_i/D_rel) ? x<=sqrt(sq(D_i/D_rel)-sq(y)) : 0)
			//f[]=f0[]*(x>=0 && y>=0 && y<=D_i/D_rel && ((y >= 0 && y <= D_i) ? x>=ell_sf_A*sqrt(sq(D_i)-sq(y)) : (y > D_i && y <= D_i/D_rel) ? x<=ell_sf_B*sqrt(sq(D_i/D_rel)-sq(y)) : 0) && x<=ell_sf_B*sqrt(sq(D_i/D_rel)-sq(y))); //quarter circles
			u.x[]= sq(D_i)>=sq(y)+sq(z) && x<D_i ? U_g : sq(D_i/D_rel)>=sq(y) + sq(z) && sq(D_i)<=sq(y)+sq(z) && x<=ell_sf_A*(sqrt(sq((D_i/D_rel-D_i)/2.)-sq(y-D_i-(D_i/D_rel-D_i)/2.))) ? f[]*U_g/U_rel : 0;
			//u.x[]=(y >= 0 && y < D_i && x<ell_sf_A*sqrt(sq(D_i)-sq(y))) ? 0 : (y >= D_i && y <= D_i/D_rel && x<=ell_sf_B*sqrt(sq(D_i/D_rel)-sq(y))) ? f[]*U_g/U_rel : 0;
		}
	}
}

//Run time logger
event logfile (i++){
	if (i==0)
			fprintf (stderr, "t dt mgp.i mgpf.p mgu.i grid->tn perf.t perf.speed\n");
	fprintf (stderr, "%g %g %d %d %d %ld %g %g\n",t,dt,mgp.i,mgpf.i,mgu.i,grid->tn,perf.t,perf.speed);
	
}

event perf_plot (t=T_END*U_rel) {
  if (getenv ("DISPLAY"))
    popen ("gnuplot -e 'set term x11 noraise title perfs' "
	   "$BASILISK/navier-stokes/perfs.plot "
	   "& read dummy; kill $!", "w");
}
//Output animation
event movie (t += 1e-2){

	//view(fov=1, tx = -0.5, width=1280);
	view (tx = -0.5);
	clear();
	draw_vof("f");
	//squares ("u.x", linear = true);
	box();
	mirror(n={0,1}){
		draw_vof("f");
		//squares ("u.y", linear = true);
		box();
	}

	save("movie.mp4");
}


//Output Paraview XML VTU file for post processing
event field_binout (t+=0.1){
	scalar lev[];
	char filename[256];
	sprintf(filename,"frac-%f.bin",t);
	FILE * fp_f=fopen(filename,"w");
	output_vtu_bin_foreach((scalar *){f}, (vector *){u}, fp_f, false);
	fclose(fp_f);
}

//Set simulation time
event end (t=0, i++, t<=T_END*U_rel){
}

//Wavelet-based adaptive grid refinement function
event adapt(i++){
	adapt_wavelet({f,u}, (double[]){femax,uemax,uemax,uemax},maxlevel);

	//unrefine (x>(domsize_coeff*0.99)*D_i);
	/*THis is a stop gap solution to prevent backflow during simulation
	by coarsening the mesh to the coarsest setting after 99% of computational domain length.
	However, backflow are usually caused by incorrect outflow boundary conditions,
	and volume might not be conserved 
	so this measure should only be used when necessary*/
}