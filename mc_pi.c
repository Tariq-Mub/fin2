/*
 * Monte Carlo pi calculation
 */

#include <gsl/gsl_rng.h>

double mc_pi(long m, gsl_rng *r);

double mc_pi(long m, gsl_rng *r)
{
   double x, y;       /* x,y coordinate, both between 0 and 1  */
   long count = 0L, n;

   for (n = 0; n < m; n++) {
      /* generate random coordinates */
      x = gsl_rng_uniform(r);
      y = gsl_rng_uniform(r);

      /* count when the point lands in the circle */
      if ((x*x + y*y) <= 1.0) {
         count += 1;
      }
   }

   return(4.0 * ((double)count)/((double)m));

} 
