#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "GeneralRelativity.h"

int main()
{
    
    //printf("If you were hooked up yo Voyager 1 since it left Neptune on the 25th of August,1989. You'd have aged %.4f seconds less than everyone on earth!\n", time_dilation(17.27*1000, 12169 * 24 * 60 * 60));
    double F = 5e14;
    double f_prime = GravitationalRedshift(F,2.1e30,5.5e6);
    printf("Emitted frequency is %.10e Hz\n", F );
    printf("Measured frequency is %.10e Hz\n", f_prime  );
    double delta = F - f_prime ;
    printf("ration change in frequency is %.10f Hz\n", delta / F );

    return 0;
}
