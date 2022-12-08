#include "saint-venant.h"

event end (i=10){
      printf("i=%d t= %g\n",i ,t);
}

int main(){
    run();
}
