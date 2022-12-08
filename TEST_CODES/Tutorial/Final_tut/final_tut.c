#include "saint-venant.h"

#define LEVEL 8

//this is a macro called LEVEL, which allows us to run the same code while varying the resolution, so that we dont need to compile and run multiple times

event init (t=0){
    foreach()
        h[]=0.1 + 1.*exp(-200.*(x*x +y*y));
}

event graphs (i++){
    stats s = statsf(h);
    fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

event images (t += 4./300.){ //setting time intervals
    output_ppm(h, linear=true); //enable bilinear interpolation

    scalar l[];
    foreach()
        l[] = LEVEL;
    static FILE * fp = fopen ("grid.ppm","w");
    output_ppm(l,fp,min=0,max=LEVEL);
}

event end (t=4){
    printf("i = %d t = %g\n", i, t);
}

event adapt (i++){
    adapt_wavelet ({h}, (double[]){4e-3}, maxlevel = LEVEL);
}

int main(){
    origin(-0.5,-0.5);
    init_grid (1 << LEVEL); 
//varying grid resolution according to the current LEVEL
    run();
}
