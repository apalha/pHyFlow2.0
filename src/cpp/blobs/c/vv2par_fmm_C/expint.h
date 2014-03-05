/* expint.h */

/* S. Engblom 2007-05-19 */

#ifndef __expint_h
#define __expint_h

#include <math.h>

/* Direct inclusion of source due to error message "Error: External
   calls are not supported" from nvcc. */

/*------------------------------------------------------------------------*/
#ifdef __CUDACC__
__host__ __device__ inline double expint(double x)
#else
inline double expint(double x)
#endif
/* Evaluates the exponential integral E_1(x) for non-negative real
   arguments x. The minmax approximations have been determined using
   Maple 9.5. */
{
  if (x <= 2.0) {
    /* For small x, E_1(x) = -log(x)+R(x), where R is
       well-behaved. Here R is approximated by a rational polynomial
       of order (6,6). The formal relative accuracy of this
       approximation is about 3.1e-17 but an additional source of
       error is the evaluation itself. */
//     const static double g1[] = {
//       0.000017396708731877516410,0.00028436097144890551495,
//       0.0050571525971836556513,0.011110873238262490464,
//       0.074118578348927945299,-0.86879488231965009498,0.54916988795210092634
//     };
//     const static double g2[] = {
//       0.0000034609955425018902702,0.00011558991363390501606,
//       0.0017158980920046998877,0.013012642421323580031,
//       0.035695607287898557846,-0.14313034833137060132,-0.95141196149931862057
//     };
//     double p1 = g1[0],p2 = g2[0];
// 
//     for (int i = 1; i < sizeof(g1)/sizeof(g1[0]); i++) {
//       p1 = x*p1+g1[i];
//       p2 = x*p2+g2[i];
//     }
    
    double p1 = 0.000017396708731877516410,p2 = 0.0000034609955425018902702;
    p1 = x*p1+0.00028436097144890551495;
    p2 = x*p2+0.00011558991363390501606;
    p1 = x*p1+0.0050571525971836556513;
    p2 = x*p2+0.0017158980920046998877;
    p1 = x*p1+0.011110873238262490464;
    p2 = x*p2+0.013012642421323580031;
    p1 = x*p1+0.074118578348927945299;
    p2 = x*p2+0.035695607287898557846;
    p1 = x*p1-0.86879488231965009498;
    p2 = x*p2-0.14313034833137060132;
    p1 = x*p1+0.54916988795210092634;
    p2 = x*p2-0.95141196149931862057;
    return p1/p2-log(x);
  }
  else {
    /* For large x, E_1(x) = exp(-x)/x*Q(1/x), where Q is
       well-behaved. Here Q is approximated by a rational polynomial
       of order (8,8) yielding an absolute error (in Q) of about
       3.2e-16. Again, the evaluation itself induces another (small)
       error. */
//     const static double g1[] = {
//       0.026122682555578989999,1.9571952984220793991,11.232557465288610841,
//       19.726465904375566681,14.513335769852306891,5.0791791313012124334,
//       0.87721505376439371882,0.071134153858027451920,0.0021326257031014296815
//     };
//     const static double g2[] = {
//       0.69562331955631940293,7.8125427293684917823,24.346870262625212927,
//       31.018938910624158297,18.900978193409730398,5.8916579081367280645,
//       0.94621658192026229210,0.073266779561126915534,0.0021326257031014303568
//     };
//     const double y = 1.0/x;
//     double p1 = g1[0],p2 = g2[0];
// 
//     for (int i = 1; i < sizeof(g1)/sizeof(g1[0]); i++) {
//       p1 = y*p1+g1[i];
//       p2 = y*p2+g2[i];
//     }
    const double y = 1.0/x;
    double p1 = 0.026122682555578989999,p2 = 0.69562331955631940293;
    p1 = y*p1+1.9571952984220793991;
    p2 = y*p2+7.8125427293684917823;
    p1 = y*p1+11.232557465288610841;
    p2 = y*p2+24.346870262625212927;
    p1 = y*p1+19.726465904375566681;
    p2 = y*p2+31.018938910624158297;
    p1 = y*p1+14.513335769852306891;
    p2 = y*p2+18.900978193409730398;
    p1 = y*p1+5.0791791313012124334;
    p2 = y*p2+5.8916579081367280645;
    p1 = y*p1+ 0.87721505376439371882;
    p2 = y*p2+0.94621658192026229210;
    p1 = y*p1+0.071134153858027451920;
    p2 = y*p2+0.073266779561126915534;
    p1 = y*p1+0.0021326257031014296815;
    p2 = y*p2+0.0021326257031014303568;
    return p1/p2*exp(-x)*y;
  }
}
/*------------------------------------------------------------------------*/
#endif /* __expint_h */
