#include "saint-venant.h"

event init (i=0){
    foreach()
        h[]=0.1 + 1.*exp(-200.*(x*x +y*y));
}

event graphs (i++){
    stats s = statsf(h);
    fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

event images (i++){
    output_ppm(h);
}
event end (i=300){
    printf("i = %d t = %g\n", i, t);
}

int main(){
    origin(-0.5,-0.5);
    run();
}

