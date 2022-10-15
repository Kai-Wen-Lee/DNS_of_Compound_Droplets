#include "saint-venant.h"
//deleted cartesian since we want to demonstrate adaptive grid refinement


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
}
event end (t=4){
    printf("i = %d t = %g\n", i, t);
}

event adapt (i++){
    adapt_wavelet ({h}, (double[]){4e-3}, maxlevel = 8);
}

int main(){
    origin(-0.5,-0.5);
    init_grid (256); //increasing resolution
    run();
}
