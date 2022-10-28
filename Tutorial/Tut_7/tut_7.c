// Gaussian Diffusion Tutorial

#include "diffusion.h"
#include "run.h"

scalar s[];

int main(){
	DT = 0.05; //time stepping parameter
	run();
}

event init (t=0){
	foreach()
		s[] = exp(-(sq(x - 0.5) + sq(y - 0.5))*10.);
}
//The foreach() loop iterates over all cells and sets the cell-centered coordinates (x,y) in the background.

event mov (t += 0.1){
	output_ppm (s, file = "s.mp4", n = 256, min = -1, max = 1);
	scalar lev[];
	foreach()
		lev[] = level;
	output_ppm (lev, file = "level.mp4", n = 126, max =6);
}
//n = resolution, min, mix = colour bar
event diff (i++){
	dt = dtnext (DT); //check if DT is too small
	// const face vector kap[] = {0.01, 0.01};
	face vector kap[]; //defining a variable diffusivity field
	foreach_face(x)
			kap.x[] = x/100.;
	foreach_face(y)
			kap.y[] = 0.;
	boundary ((scalar*) {kap}); //calling boundary() function to ensure proper deifinition near resolutions boundaries
	diffusion (s, dt, kap);
}

/* In the following event the time integration is performed. Again, inspired by a relevant example, we tell the time-loop to set the actual timestep (dt), based on our maximum value (DT). The documentation for the diffusion solver reveals the proper sequence of the arguments to the duffusion() function (albeit via a structure). A constant diffusivity field (kap) is declared and initialized. The values are defined on cell faces. The faces have distinct directions and hence we need to set the corresponding vector components in both dimensions.*/

event lot (i += 5){
	static FILE * fp = fopen ("data", "w");
	fprintf (fp, "%g %d %.8g\n", t, i, statsf(s).sum);
}

//Finally, a data file is writen to check if the scalar field s is conserved.

event adapt (i++)
	adapt_wavelet ({s}, (double[]){0.01}, 15);

event stop (t = 10);

/* The time loop stops when it “sees” no further events. Therefore, we request it to continue to go to the stop event at t = 10. This event could have many other names.*/