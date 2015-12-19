/* 
 * Monte Carlo calculation Pi
 */
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

/* mc_pi.c */
double mc_pi (long m, gsl_rng * r);
double g (double *t, size_t dim, void *params);

//#define POINTS 524288           /* 2^19, initial number of random points to generate */
//#define NEXP 64                 /* number of experiments for each value of points */
//#define M 11                    /* number of different points values */

int main (void)
{
    gsl_rng *r;
    unsigned long seed = 1UL;
    size_t dim = 6;
    int calls = 1000000;
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
        printf ("%10.8f   %10.8f    %10.8f\n", dist, energy, -2./pow(dist, 3.));
        fflush (stdout);
        dist += dstep;
    }

    gsl_rng_free (r);

    return 0;
}
