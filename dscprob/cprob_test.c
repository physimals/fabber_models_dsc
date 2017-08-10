#include "cprob.h"

#include <stdio.h>

#define VAL 3.533
#define VAL2 0.731

/**
 * This is a really crude random set of sanity
 * checks. Have had issues with inconsistent
 * results in these functions due to compile
 * errors and existance of non-standard math
 * library functions, so just wrote this
 * as a simple way to see if the values are
 * coming out the same on different platforms
 */
int main(int argc, char **argv) 
{
    printf("gamma(%f) = %e\n", VAL, gamma(VAL));
    printf("lgam(%f) = %e\n", VAL, lgam(VAL));
    printf("expm1(%f) = %e\n", VAL, expm1(VAL));
    printf("log1p(%f) = %e\n", VAL, log1p(VAL));
    printf("erf(%f) = %e\n", VAL, erf(VAL));
    printf("erfc(%f) = %e\n", VAL, erfc(VAL));
    printf("igam(%f, %f) = %e\n", VAL, VAL2, igam(VAL, VAL2));
    printf("igami(%f, %f) = %e\n", VAL, VAL2, igami(VAL, VAL2));
}
