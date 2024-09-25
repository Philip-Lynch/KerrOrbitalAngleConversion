#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <math.h>

#include "constants_of_motion.h"
#include "angle_conversions.h"

// Structure to hold the parameters
typedef struct {
    double a;
    double p;
    double e;
    double x;
    double Phi_r;
    double Phi_theta;
} Params;

// Define the system of equations
int system_of_equations(const gsl_vector *q, void *params, gsl_vector *f) {
    double qr = gsl_vector_get(q, 0); 
    double qz = gsl_vector_get(q, 1);

    Params *point = (Params *)params;      // Cast void pointer to Params pointer
    double a = point->a;
    double p = point->p;
    double e = point->e;
    double x = point->x;
    double Phi_r = point->Phi_r;
    double Phi_theta = point->Phi_theta;

    double f0 = Phi_r - RadialMinoPhaseToBoyerLindquistPhase(a,p,e,x,qr,qz);
    double f1 = Phi_theta - PolarMinoPhaseToBoyerLindquistPhase(a,p,e,x,qr,qz);

    gsl_vector_set(f, 0, f0);
    gsl_vector_set(f, 1, f1);

    return GSL_SUCCESS;
}

// Function to solve the system of equations
int RootFindingMinoPhases(double a, double p,double e, double x_inc, double Phi_r, double Phi_theta, double *qr, double *qz) {
    
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *s;
    int status;
    size_t iter = 0;
    double tolerance = 1e-10;

    // Number of equations and unknowns
    const size_t n = 2;

    // Initial guess for the roots 
    gsl_vector *q = gsl_vector_alloc(n);
    gsl_vector_set(q, 0, Phi_r);  // Initial guess for qr
    gsl_vector_set(q, 1, Phi_theta);  // Initial guess for qz

    // Parameters a and b
    Params params = {a,p,e,x_inc,Phi_r,Phi_theta};

    // Define the system of equations
    gsl_multiroot_function f = {&system_of_equations, n, &params};

    // Allocate the solver
    s = gsl_multiroot_fsolver_alloc(T, n);
    gsl_multiroot_fsolver_set(s, &f, q);

    // Iterate until the solver converges or fails
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        if (status)   // Check if solver is stuck
            break;

        status = gsl_multiroot_test_residual(s->f, tolerance);  // Test for convergence

    } while (status == GSL_CONTINUE && iter < 1000);

    if (status == GSL_SUCCESS) {
        // Store results in the output parameters
        *qr = gsl_vector_get(s->x, 0);
        *qz = gsl_vector_get(s->x, 1);
    } else {
        printf("Failed to converge\n");
        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(q); 
        return status;
    }

    // Clean up
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(q);

    return GSL_SUCCESS;
}
