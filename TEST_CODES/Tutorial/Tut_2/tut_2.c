#include "saint-venant.h"

event init (t=0){ //set up initial conditions
  foreach()
    h[]=0.1 + 1.*exp(-200.*(x*x +y*y));
}

event images (i++){ //print simple images
    output_ppm(h);
}

event end (i=300){
  printf("i = %d t = %g\n", i, t);
}

int main(){
  origin (-0.5,-0.5); //Setting point of origin
  run();
}
