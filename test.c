#include <stdio.h>

#include "poly.h"

#define DEGREE 5

#define NO_POINTS 11

#define LENGTH( x ) ( sizeof( x ) / sizeof( x[0] ))

double x[] = {0,  1,  2,  3,  4,  5,  6,   7,   8,   9,   10};
double y[] = {1,  6,  17, 34, 57, 86, 121, 162, 209, 262, 321};

double p1[] = {2, 2, 2, 2};
double p2[] = {20, 40};
double p3[] = {1, 0, 0};    // x^2
double x0 = 1;
double x1 = 2;

int main(int argc, char ** argv){
    uint32_t i;
    printf("\n-- FITTING ----------------------------------\n");
    poly_t * poly_points = polyfit(x, y, NO_POINTS, DEGREE);

    if(poly_points){
        polyprint(poly_points, 4, false, false);
        printf("r ^ 2 = %f\n", polyfitness(poly_points, x, y, NO_POINTS));

        for(i = 0; i < NO_POINTS; i++){
            printf("%f\t%f\t%f\n", x[i], y[i], polyval(poly_points, x[i]));
        }
    }

    printf("\n-- DERIVATION ----------------------------------\n");
    poly_t * poly1 = polycreate(p1, LENGTH( p1 ));
    if(poly1){
        polyprint(poly1, 1, false, true);
        for(i = 0; i < 4; i++){
            polyderivs(poly1);
            printf(" d%d : ", i + 1);
            polyprint(poly1, 1, false, true);
        }
    }
    printf("\n-- INTEGRATION ----------------------------------\n");
    poly_t * poly2 = polycreate(p2, LENGTH( p2 ));
    if(poly2){
        polyprint(poly2, 1, false, true);
        for(i = 0; i < 4; i++){
            polyintegs(poly2, ((i + 1) * 10));
            printf(" i%d : ", i + 1);
            polyprint(poly2, 1, false, true);
        }
    }

    printf("\n\n");
    poly_t * poly3 = polycreate(p3, LENGTH( p3 ));
    if(poly2){
        polyprint(poly3, 1, false, false);
        printf(" from %.1f to %.1f = %f\n", x0, x1, polyintegdef(poly3, x0, x1));
    }

    polyfree(poly_points);
    polyfree(poly1);
    polyfree(poly2);
    polyfree(poly3);

    return 0;
}
