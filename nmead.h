/**** NELDER MEAD OPTIMISATION ****/
// Adapted from https://github.com/matteotiziano/nelder-mead

/*****************************************************************************
 * Some documentation:
 *
 * To use this Nelder-Mead implementation, the use must write a cost function
 * and pass its handle to the nelder_mead() function, whose prototype is:
 *
 *   void nelder_mead(
 *                     double *,
 *                     int,
 *                     optimset_t,
 *                     point_t *,
 *                     double (*cost_fun)(int, const double *, void *),
 *                     void *
 *                   );
 *
 * The arguments are, in order:
 *
 *   double      *x0        : the initial guess of n parameter values
 *   int          n         : the number of parameters
 *   optimset_t   optimset  : a special struct, defined in this header file,
 *                            containing values specific to the Nelder-Mead
 *                            algorithm (see below)
 *   point_t     *solution  : where the final solution is stored
 *   double (*cost_fun)(...): the function pointer to the cost function
 *   void     *cost_fun_args: this pointer is passed to the cost function as
 *                            the third argument
 *
 * The form of the cost function is
 *
 *   double cost_fun(
 *                    int,
 *                    const double *,
 *                    void *
 *                  );
 *
 * Its arguments are, in order:
 *
 *   int n          : the number of fitting parameters
 *   const double * : a pointer to the n parameters
 *   void *         : a pointer to the data to be fitted
 *
 ****************************************************************************/


#ifndef NMEAD_H
#define NMEAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// define the point in the simplex (x) and its functional value (fx)
typedef struct nm_point_t
{
    double *x;
    double fx;
} nm_point;

// define optimization settings
typedef struct nm_optimset_t
{
    double tolx;
    double tolf;
    int max_iter;
    int max_eval;
} nm_optimset;

#define NM_TOL_X    0.001 // tolerance on the simplex solutions coordinates
#define NM_TOL_F    0.001 // tolerance on the function value
#define NM_MAX_ITER 1000  // maximum number of allowed iterations
#define NM_MAX_EVAL 1000  // maximum number of allowed function evaluations

#define NM_RHO      1
#define NM_CHI      2
#define NM_GAMMA    0.5
#define NM_SIGMA    0.5

#define SQUARE(x) ((x)*(x))

// Nelder-Mead algorithm and template cost function 
void nelder_mead( double *,
                  int,
                  nm_optimset,
                  nm_point *,
                  double (*cost_fun)(int, const double *, void *),
                  void *);

#endif
