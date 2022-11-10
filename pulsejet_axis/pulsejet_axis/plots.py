# -*- coding: utf-8 -*-
# #include "axi.h"
# #include "navier-stokes/centered.h"
# #include "two-phase.h"
# #include "tension.h"
# #include "tag.h"
# #include "view.h"
# 
# //Define radius of jet, init jet len, Re and surface tention coeff
# 
# #define radius_outer 12.56/100.
# #define radius_inner 9.02/100.
# #define length 0.0005
# #define Re 500
# #define SIGMA 3e-5
# #define Fr 0.88
# 
# int maxlevel=10;	//max level of refinement =10
# double uemax=0.1;	//error threshold of velocity is 0.1
# 
# /* -------------------imposing boundary conditions------------------------------*/
# scalar f0[];												 //set an aux volume fraction field (1 is inside cylinder, 0 is outside cylinder)
# u.n[left] = dirichlet(f0[]);		//applying oscilating inflow velocity on lhs
# u.t[left] = dirichlet(0); 									//applying BC of zero tangential velocity at inflow
# 
# #if dimension>2
# u.r[left] = dirichlet(0); 		//applying bc of zero radial velocity at left boundary
# #endif
# 
# 
# p[left] = neumann(0); 			//applying bc of zero pressure gradient at left boundary
# f[left] = f0[];
# 
# u.n[right] = neumann(0); 			//applying bc of zero velocity gradient at right boundary
# p[right] = dirichlet(0);			//applying bc of zero pressure at right boundary (?) 
# /* ----------------------main program---------------------------------*/
# int main(int argc, char * argv[]){
# 	if (argc>1)
# 		maxlevel = atoi (argv[1]);	//optional cl arguments: maximum level
# 	if (argc >2)
# 		uemax = atof (argv[2]);		//optional cl arguments: error threshold
# 
# 	init_grid(64); //discretise initial domain with 64^3 grid points
# 	size(10.*radius_outer);
# 	rho1=1000.,rho2=1.225;
# 	mu1=2.*radius_outer/Re*rho1;
# 	mu2=2.*radius_inner/Re*rho2;
# 	f.sigma=SIGMA;
# 
# 	run();
# 
# }
# 
# /* -------------------setting initial conditions -------------------------------*/
# event acceleration (i++) {
#   face vector av = a;
#   foreach_face(x)
#   av.x[] += 0.92;
# }
# 
# event init (t=0){
# 	if(!restore (file="restart")){
# 
# 		/*use a static refinement down to maxlevel in a cyl. 1.2 times longer than 
# 		the init jet and twice the radies*/
# 		refine (x<20*length && y<2.*radius_outer && level <maxlevel);
# 
# 
# 		/* init the aux vol fraction field held for a cyl of constant radius
# 		double outer_circle = y - radius_outer;
# 		double inner_circle = y - radius_inner;
# 		double annulus = difference (outer_circle, inner_circle);*/
# 		fraction (f0, difference((radius_outer - y), (radius_inner - y)));
# 		f0.refine=f0.prolongation=fraction_refine;
# 		restriction ({f0});		//for BC on levels
# 
# 		/* use this to define the init jet and its velocity*/
# 		foreach(){
# 			f[]=f0[]*(x<length);
# 			u.x[]=f[];
# 		}
# 	}
# }
# 
# /* -------------------Output -------------------*/
# 
# event logfile (i++){
# 	if (i==0)
# 			fprintf (stderr, "t dt mgp.i mgpf.p mgu.i grid->tn perf.t perf.speed\n");
# 	fprintf (stderr, "%g %g %d %d %d %ld %g %g\n",t,dt,mgp.i,mgpf.i,mgu.i,grid->tn,perf.t,perf.speed);
# }
# 
# //generating animation using basilisk view
# event movie (t += 1e-2){
# 
# 	view (tx = -5.*radius_outer);
# 	clear();
# 	draw_vof ("f");
# 	box();
# 
# 	save("movie.mp4");
# }
# 
# event snapshot (t=0.1;t+=0.1;t<=3.8){
# 	char name[80];
# 	sprintf (name, "snapshot-%g", t);
# 	scalar pid[];
# 	foreach()
# 		pid[]=fmod(pid()*(npe()+37),npe());
# 	dump (name);
# }
# 
# event adapt(i++){
# 	adapt_wavelet({f,u}, (double[]){0.01,uemax,uemax,uemax},maxlevel);
# }
