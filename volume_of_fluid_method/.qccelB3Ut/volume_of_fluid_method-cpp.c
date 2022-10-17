@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "volume_of_fluid_method.c"
//Two phase interfacial flows

#include "vof.h"

/* volume fraction in fluid 1 is f=1 and f=0 in fluid 2
densities and dynamic vsicoucity for fluid 1 and 2 are rho1,mu1,rho2,mu2 */

scalar f[], * interfaces = {f};
double rho1=1., mu1=0., rho2=1., mu2=0.;

/* aux fields are necessary to define the variable specific volume alpha = 1/rho
as well as the cell-centered density */

face vector alphav[];
scalar rhov[];

event defaults (i=0){
	alpha=alphav;
	rho=rhov;

	/* if viscosity is non-zero, we need to allocate the face-centered viscosity field. */

	if (mu1 || mu2)
		mu = new face vector;
		
	/* Add interface to the default display */
	
	display ("draw_vof (c='f');");
}

/* defining density and viscocity
defined using arithmetic averages by default */

#ifndef rho
#define rho(f) (clamp(f,0.,1.)*(rho1-rho2)+rho2)
#endif

#ifndef mu
#define mu(f) (clamp(f,0.,1.)*(mu1-mu2)+mu2)
#endif

/* "smearing" of the density/viscosity jump */
#ifdef FILTERED
scalar sf[];
#else
#define sf f
#endif

event tracer_advection (i++){
#ifndef sf
#if dimension <= 2
	foreach()
		sf[]=(4.*f[]+2.*(f[0,1]+f[0,-1]+f[-1,0])+f[-1,-1]+f[1,-1]+f[1,1]+f[-1,1])/16.;
#else
	foreach()
		sf[]=(8.*f[]+4.*(f[-1]+f[1]+f[0,1]+f[0,-1]+f[0,0,1]+f[0,0,-1])+
					2.*(f[-1,1]+f[-1,0,1]+f[-1,0,-1]+f[-1,-1]+
						f[0,1,1]+f[0,1,-1]+f[0,-1,1]+f[0,-1,-1]+
						f[1,1]+f[1,0,1]+f[1,-1]+f[1,0,-1])+
					f[1,-1,1]+f[-1,1,1]+f[-1,1,-1]+f[1,1,1]+
					f[1,1,-1]+f[-1,-1,-1]+f[1,-1,-1]+f[-1,-1,1])/64.;
#endif
#endif
	
#if TREE
	sf.prolongation = refine;
	sf.dirty=true;
#endif
}

event properties (i++)
{
	foreach_face(){
		double ff = (sf[]+sf[-1])/2.;
		alphav.x[] = fm.x[]/rho(ff);
		if (mu1||mu2){
			face vector muv=mu;
			muv.x[]=fm.x[]*mu(ff);
			}
		}
	foreach()
		rhov[]=cm[]*rho(sf[]);
#if TREE
	sf.prolongation=fraction_refine;
	sf.dirty = true;
#endif
}

int main(){
	run();
}

#endif
