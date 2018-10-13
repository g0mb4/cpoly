#include <stdio.h>

#include "poly.h"

#define DEGREE 5

#define NO_POINTS 11

double x[] = {0,  1,  2,  3,  4,  5,  6,   7,   8,   9,   10};
double y[] = {1,  6,  17, 34, 57, 86, 121, 162, 209, 262, 321};

int main(int argc, char ** argv){
    uint32_t i;
    poly_t * poly = polyfit(x, y, NO_POINTS, DEGREE);

    if(poly){
        polyprint(poly, 4, false);
        printf("r ** 2 = %f\n", polyfitness(poly, x, y, NO_POINTS));

        for(i = 0; i < NO_POINTS; i++){
            printf("%f\t%f\t%f\n", x[i], y[i], polyval(poly, x[i]));
        }
    }

    polyfree(poly);

    return 0;
}
