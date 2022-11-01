#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

//Define radius of jet, init jet len, Re and surface tention coeff

#define radius 1./4.
#define length 0.025
#define Re 5800
#define SIGMA 3e-5

int maxlevel=10;	//max level of refinement =10
double uemax=0.1;	//error threshold of velocity is 0.1

/* -------------------imposing boundary conditions------------------------------*/
scalar f0[];												 //set an aux volume fraction field (1 is inside cylinder, 0 is outside cylinder)
u.n[left] = dirichlet(f0[]*(1.+0.5*sin(10.*2.*pi*t)));		//applying oscilating inflow velocity on lhs
u.t[left] = dirichlet(0); 									//applying BC of zero tangential velocity at inflow

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
	size(3.);

	rho1=1.,rho2=1./27.84;
	mu1=2.*radius/Re*rho1;
	mu2=2.*radius/Re*rho2;
	f.sigma=SIGMA;

	run();

}

/* -------------------setting initial conditions -------------------------------*/
event init (t=0){
	if(!restore (file="restart")){

		/*use a static refinement down to maxlevel in a cyl. 1.2 times longer than 
		the init jet and twice the radies*/
		refine (x<1.2*length && y<2.*radius && level <maxlevel);


		/* init the aux vol fraction field held for a cyl of constant radius*/
		fraction (f0,radius-y);
		f0.refine=f0.prolongation=fraction_refine;
		restriction ({f0});		//for BC on levels

		/* use this to define the init jet and its velocity*/
		foreach(){
			f[]=f0[]*(x<length);
			u.x[]=f[];
		}
	}
}

/* -------------------Output -------------------*/

event logfile (i++){
	if (i==0)
			fprintf (stderr, "t dt mgp.i mgpf.p mgu.i grid->tn perf.t perf.speed\n");
	fprintf (stderr, "%g %g %d %d %d %ld %g %g\n",t,dt,mgp.i,mgpf.i,mgu.i,grid->tn,perf.t,perf.speed);
}

//generating animation using basilisk view
event movie (t += 1e-2){
#if dimension ==2
	scalar omega[];
	vorticity (u,omega);
	view (tx =-0.5);
	clear();
	draw_vof ("f");
	squares ("omega", linear = true, spread =10);
	box();
#else
	scalar pid[];
	foreach()
			pid[]=fmod(pid()*(npe()+37),npe());
	view (camera ="iso", fov=14.5,tx=-0.418,ty=0.288,width=1600,height=1200);
	clear();
	draw_vof("f");
#endif
	save("movie.mp4");
}

event snapshot (t=0.1;t+=0.1;t<=3.8){
	char name[80];
	sprintf (name, "snapshot-%g", t);
	scalar pid[];
	foreach()
		pid[]=fmod(pid()*(npe()+37),npe());
	dump (name);
}

event adapt(i++){
	adapt_wavelet({f,u}, (double[]){0.01,uemax,uemax,uemax},maxlevel);
}
