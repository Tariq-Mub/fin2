/* 
 * Monte Carlo calculation Pi
 */
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

double g (double *t, size_t dim, void *params);

int main (void)
{
    gsl_rng *r;
    unsigned long seed = 1UL;
    size_t dim = 6;
    int calls = 10000000;
    double t[6];
    int i, j, l;
    int np = 20;
    double dmin = 1.001, dmax = 4., dist;
    double dstep = (dmax - dmin) / (np - 1);
    double energy;

    /* allocate random number generator */
    r = gsl_rng_alloc (gsl_rng_taus2);

    /* seed the random number generator */
    gsl_rng_set (r, seed);

    dist = dmin;
    for (i = 0; i < np; i++)
    {
        double sums = 0;

        for (l = 0; l < calls; l++)
        {
            for (j = 0; j < (int) dim; j++)
            {
                t[j] = gsl_rng_uniform (r);
            }
            sums += g (t, dim, &dist);
        }
        energy = sums / calls;
        printf ("%10.8f   %10.8f    %10.8f\n", dist, energy, -2. / pow (dist,
                3.));
        fflush (stdout);
        dist += dstep;
    }

    double xl[] = { 0., 0., 0., 0., 0., 0. };
    double xu[] = { 1., 1., 1., 1., 1., 1. };
    double res, err;

    gsl_monte_function G = { &g, dim, &dist };

    gsl_monte_vegas_state *sv = gsl_monte_vegas_alloc (dim);

    gsl_monte_vegas_init (sv);

    printf ("# Dist    Result ErrEst Dipol\n");
    dist = dmin;
    for (i = 0; i < np; i++)
    {

        gsl_monte_vegas_integrate (&G, xl, xu, dim, (size_t) calls / 5, r, sv, &res,
            &err);

        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, dim, (size_t) calls, r, sv, &res,
                &err);
        }
        while (fabs (gsl_monte_vegas_chisq (sv) - 1.0) > 0.2);
        printf ("%10.8f   %10.8f    %10.8f    %10.8f\n", dist, res, err,
            -2. / pow (dist, 3.));
        fflush (stdout);
        dist += dstep;
    }

    gsl_monte_vegas_free (sv);

    gsl_rng_free (r);

    return 0;
}
