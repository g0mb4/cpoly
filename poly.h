#ifndef __POLY_H__
#define __POLY_H__

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include <gsl/gsl_multifit.h>

typedef struct S_POLY {
    uint32_t degree;
    double * coeffs;
} poly_t;

poly_t * polycreate(double * coeffs, uint32_t size);
void polyfree(poly_t * poly);
void polycopy(poly_t * dest, poly_t * src);

poly_t * polyfit(double * px, double * py, size_t size, uint32_t degree);
double polyval(poly_t * poly, double x);
double polyfitness(poly_t * poly, double * px, double * py, size_t size);

poly_t * polyderiv(poly_t * poly);
void polyderivs(poly_t *poly);

poly_t * polyinteg(poly_t * poly, double C);
void polyintegs(poly_t * poly, double C);
double polyintegdef(poly_t * poly, double x0, double x1);

void polyprint(poly_t * poly, uint32_t precision, bool full, bool degree);

#endif
