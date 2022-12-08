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
#line 1 "tut_2.c"
#include "saint-venant.h"

event init (t=0){ //set up initial conditions
  foreach(){
    h[]=0.1 + 1.exp*(-200.*(x*x +y*y));}
}

event end (i=10){
  printf("i = %d t = %g\n", i, t);
}

int main(){
  origin (-0.5,-0.5); //Setting point of origin
  run();
}

#endif
