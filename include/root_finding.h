#include <gsl/gsl_vector.h>

int system_of_equations(const gsl_vector *q, void *params, gsl_vector *f);

int RootFindingMinoPhases(double a, double p,double e, double x, double Phi_r, double Phi_theta, double *qr, double *qz);