#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

double g (double *t, size_t dim, void *params);

int main (void)
{

    size_t dim = 4;
    double xl[] = { -1., -1., -1., -1. };
    double xu[] = { 1., 1., 1., 1. };

    gsl_rng *r = gsl_rng_alloc (gsl_rng_taus2);
    unsigned long seed = 1UL;

    gsl_rng_set (r, seed);

    size_t calls = 10000000;

    gsl_monte_function G = { &g, dim, NULL };

    gsl_monte_vegas_state *sv = gsl_monte_vegas_alloc (dim);

    gsl_monte_vegas_init (sv);

    gsl_monte_vegas_integrate (&G, xl, xu, dim, calls / 5, r, sv, &res, &err);

    printf ("# Stat Result ErrEst Exact AbsErr\n");
    printf (" %.2f % .6f % .6f %.6f %.6f\n", gsl_monte_vegas_chisq (sv), res,
        err, exact, fabs (res - exact));
    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, dim, calls, r, sv, &res, &err);
        printf (" %.2f % .6f % .6f %.6f %.6f\n", gsl_monte_vegas_chisq (sv),
            res, err, exact, fabs (res - exact));
        fflush (stdout);
    }
    while (fabs (gsl_monte_vegas_chisq (sv) - 1.0) > 0.2);
    printf (" % .6f % .6f %.6f %.6f\n", res, err, exact, fabs (res - exact));

    gsl_monte_vegas_free (sv);

    return 0;
}
