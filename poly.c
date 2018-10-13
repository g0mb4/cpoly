#include "poly.h"

poly_t * polycreate(double * coeffs, uint32_t size){
    uint32_t i;
    poly_t * poly = (poly_t *)malloc(sizeof(poly_t));
    if(poly == NULL)
        return NULL;

    if(size > 0){
        poly->degree = size - 1;

        poly->coeffs = (double *)malloc(sizeof(double) * (size));
        if(poly->coeffs == NULL){
            free(poly);
            return NULL;
        }

        for(i = 0; i < size; i++){
            poly->coeffs[i] = coeffs[i];
        }

        return poly;
    } else {
        return NULL;
    }
}

void polyfree(poly_t * poly){
    if(poly){
        if(poly->coeffs){
            free(poly->coeffs);
        }
        free(poly);
    }
}

/*
 source : http://rosettacode.org/wiki/Polynomial_regression#C
 */
poly_t * polyfit(double * px, double * py, size_t size, uint32_t degree){
    gsl_multifit_linear_workspace * ws;
    gsl_matrix * cov, * X;
    gsl_vector * y, * c;
    double chisq;

    uint32_t i, j;

    poly_t * poly = (poly_t *)malloc(sizeof(poly_t));
    if(poly == NULL)
        goto error;

    poly->degree = degree;
    degree += 1;    // second deg: a x ** 2 + b x ** 1 + c => 3 coeffs

    poly->coeffs = (double *)malloc(sizeof(double) * (degree));
    if(poly->coeffs == NULL){
        goto error;
    }

    X = gsl_matrix_alloc(size, degree);
    if(X == NULL){
        goto error;
    }

    y = gsl_vector_alloc(size);
    if(y == NULL){
        goto error;
    }

    c = gsl_vector_alloc(degree);
    if(c == NULL){
        goto error;
    }

    cov = gsl_matrix_alloc(degree, degree);
    if(cov == NULL){
        goto error;
    }

    for(i = 0; i < size; i++) {
        for(j = 0; j < degree; j++) {
            gsl_matrix_set(X, i, j, pow(px[i], j));
        }
        gsl_vector_set(y, i, py[i]);
    }

    ws = gsl_multifit_linear_alloc(size, degree);
    gsl_multifit_linear(X, y, c, cov, &chisq, ws);

    for(i = 0; i < degree; i++) {
        poly->coeffs[degree - i - 1] = gsl_vector_get(c, i);
    }

    goto memfree;

error:
    polyfree(poly);

memfree:
    if(ws)
        gsl_multifit_linear_free(ws);

    if(X)
        gsl_matrix_free(X);

    if(cov)
        gsl_matrix_free(cov);

    if(y)
        gsl_vector_free(y);

    if(c)
        gsl_vector_free(c);

    return poly;
}

double polyval(poly_t * poly, double x){
    double sum = 0;
    if(poly){
        uint32_t i;
        for(i = 0; i < poly->degree + 1; i++){
            sum += poly->coeffs[i] * pow(x, poly->degree - i);
        }
    }
    return sum;
}

double polyfitness(poly_t * poly, double * px, double * py, size_t size){
    if(poly == NULL){
        return 0;
    }

    uint32_t i;
    double ybar = 0;
    for(i = 0; i < size; i++){
        ybar += py[i];
    }
    ybar /= (double)size;

    double ssreg = 0;
    for(i = 0; i < size; i++){
        ssreg += pow(polyval(poly, px[i]) - ybar, 2);
    }

    double sstot = 0;
    for(i = 0; i < size; i++){
        sstot += pow(py[i] - ybar, 2);
    }

    return ssreg / sstot;   // R^2
}

poly_t * polyderiv(poly_t * poly){
    poly_t * ret = (poly_t *)malloc(sizeof(poly_t));
    if(ret == NULL)
        return NULL;

    /* TODO */

    return ret;
}

poly_t * polyinteg(poly_t * poly, double C){
    poly_t * ret = (poly_t *)malloc(sizeof(poly_t));
    if(ret == NULL)
        return NULL;

    /* TODO */

    return ret;
}

void polyprint(poly_t * poly, uint32_t precision, bool full){
    if(poly){
        uint32_t i;
        char ch;
        char fmt[32];
        char fmt_last[32];

        precision = precision > 15 ? 15 : precision;

        sprintf(fmt, "%%c %%.%df x ^ %%d ", precision);
        sprintf(fmt_last, "%%c %%.%df ", precision);

        for(i = 0; i < poly->degree + 1; i++){
            ch = poly->coeffs[i] < 0 ? '-' : '+';

            if(full){
                printf(fmt, ch, fabs(poly->coeffs[i]), poly->degree - i);
            } else {
                if(fabs(poly->coeffs[i]) < 1.0 / pow(10, precision)){
                    continue;
                } else if(i == poly->degree){
                    printf(fmt_last, ch, fabs(poly->coeffs[i]));
                } else {
                    printf(fmt, ch, fabs(poly->coeffs[i]), poly->degree - i);
                }
            }
        }
        printf("\n");
    } else {
        printf("invalid polynome\n");
    }
}
