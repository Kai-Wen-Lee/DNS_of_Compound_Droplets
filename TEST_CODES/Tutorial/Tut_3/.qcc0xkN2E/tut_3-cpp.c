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
#line 1 "tut_3.c"
#include "saint-venant.h"

event init (t=0){
    foreach()
        h[]=0.1 + 1.*exp(-200.*(x*x +y*y));
}

event graphs (i++){
    stats s = statsf(h);
    forintf (stderr, "%g %g %g\n", t, s.min, s.max");
}

event images (i++){
    output_ppm(h);
}
event end (t=300){
    printf("i = %d t = %g\n", i, t);
}

int main(){
    origin(-0.5,-0.5);
    run();
}


#endif
