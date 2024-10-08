
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <math.h>
#include<float.h>
#include "common_definitions.h"

#define N 30



// Elliptic integrals use Mathematica's definitions

double EllipticK(double k){
	return gsl_sf_ellint_Kcomp(sqrt(k), GSL_PREC_DOUBLE);
}

double EllipticF(double phi, double k){
	return gsl_sf_ellint_F(phi, sqrt(k), GSL_PREC_DOUBLE) ;
}

double EllipticE(double k){
	return gsl_sf_ellint_Ecomp(sqrt(k), GSL_PREC_DOUBLE);
}

double EllipticEIncomp(double phi, double k){
	return gsl_sf_ellint_E(phi, sqrt(k), GSL_PREC_DOUBLE) ;
}

double EllipticPi(double n, double k){
	return gsl_sf_ellint_Pcomp(sqrt(k), -n, GSL_PREC_DOUBLE);
}

double EllipticPiIncomp(double n, double phi, double k){
	return gsl_sf_ellint_P(phi, sqrt(k), -n, GSL_PREC_DOUBLE);
}

/*
double Jacobi_dn(double u, void * params)
{
  int m = *(int *) params;
  double f = gsl_pow_int(x, m) + 1.0;
  return f;
}

double Jacobi_Amp(double u, double m){
double sn, cn, dn, result, 


return result
}

*/
double Jacobi_am(double u, char arg,  double x)
{
   long double a[N+1];
   long double g[N+1];
   long double c[N+1];
   long double two_n;
   long double phi;
   long double k;
   int n;

                        // Check special case x = 0 //
                        // i.e. k = m = alpha = 0.  //

   if ( x == 0.0 ) return u;

   switch (arg) {
      case 'a': k = sinl( fabsl((long double) x) ); break;
      case 'm': k = sqrtl( fabsl((long double) x) ); break;
      default:  k = (long double) fabs(x);
   }

                   // Check special case k = 1 //

   if ( k == 1.0 ) return 2.0 * atan( exp(u) ) - M_PI_2;

         // If k > 1, then perform a Jacobi modulus transformation. //
         // Initialize the sequence of arithmetic and geometric     //
         // means, a = 1, g = k'.                                   //

   a[0] = 1.0L;
   g[0] = sqrtl(1.0L - k * k);
   c[0] = k;
   
   // Perform the sequence of Gaussian transformations of arithmetic and //
   // geometric means of successive arithmetic and geometric means until //
   // the two means converge to a common mean (upto machine accuracy)    //
   // starting with a = 1 and g = k', which were set above.              //
   
   two_n = 1.0L; 
   for (n = 0; n < N; n++) {
      if ( fabsl(a[n] - g[n]) < (a[n] * LDBL_EPSILON) ) break;
      two_n += two_n;
      a[n+1] = 0.5L * (a[n] + g[n]);
      g[n+1] = sqrtl(a[n] * g[n]);
      c[n+1] = 0.5L * (a[n] - g[n]);
   }

         // Prepare for the inverse transformation of phi = x * cm. //

   phi = two_n * a[n] * u;

                      // Perform backward substitution //

   for (; n > 0; n--) phi = 0.5L * ( phi + asinl( c[n] * sinl(phi) / a[n]) );

   return (double) phi; 
}

// Inner and Outer Horizons

double InnerHorizon(double a){

    double M = 1.0;

    return M - Sqrt(-Power(a,2) + Power(M,2));

}

double OuterHorizon(double a){

    double M = 1.0;

    return M + Sqrt(-Power(a,2) + Power(M,2));
}

// Primary roots of motion
double RadialRoot1(double p, double e){
    return p/(1 - e);
}


double RadialRoot2(double p, double e){
    return p/(1 + e);
}

double PolarRoot1(double x){
    return Sqrt(1 - Power(x,2));
}

//Needed for E,L and Q calculaions

double Delta(double a, double r){
        return Power(a,2) - 2.*r + Power(r,2);
}

double D(double a, double r, double zm){
        return (Power(r,2) + Power(a,2)*Power(zm,2))*Delta(a,r);
}

double F(double a, double r, double zm){
        return Power(r,4) + Power(a,2)*(r*(2. + r) + Power(zm,2)*Delta(a,r));
}

double G(double a, double r){
        return 2. * a * r;
}

double H(double a, double r, double zm){
        return (-2. + r)*r + (Power(zm,2)*Delta(a,r))/(1. - 1.*Power(zm,2));
}


// Constants of Motion



double Energy(double a, double p, double e, double x){

    double kappa, epsilon, rho, eta, sigma,r1,r2,zm;

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    zm = PolarRoot1(x);

    kappa = D(a,r1,zm)*H(a,r2,zm) - H(a,r1,zm)*D(a,r2,zm);
    epsilon = D(a,r1,zm)*G(a,r2) - G(a,r1)*D(a,r2,zm);
    rho = F(a,r1,zm)*H(a,r2,zm) - H(a,r1,zm)*F(a,r2,zm);
    eta = F(a,r1,zm)*G(a,r2) - G(a,r1)*F(a,r2,zm);
    sigma = G(a,r1)*H(a,r2,zm) - H(a,r1,zm)*G(a,r2);

    return Sqrt((kappa*rho + 2.*epsilon*sigma - 2.*x*Sqrt(sigma*(-1.*eta*Power(kappa,2) 
    + epsilon*kappa*rho + Power(epsilon,2)*sigma)/Power(x,2)))/(Power(rho,2) + 4.*eta*sigma));
}

double AngularMomentum(double a, double p, double e, double x){

    double En, r1, r2, zm;
    
    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    zm = PolarRoot1(x);
    
    En = Energy(a,p,e,x);

    return (-1.*En*G(a,r1) + x*
      Sqrt((-1.*D(a,r1,zm)*H(a,r1,zm) + Power(En,2)*(Power(G(a,r1),2) 
      + F(a,r1,zm)*H(a,r1,zm)))/(Power(x,2))))/H(a,r1,zm);
}

double CarterConstant(double a, double p, double e, double x){
    
    double En, L,zm;

    zm = PolarRoot1(x);
    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);

    return Power(zm,2)*(Power(a,2)*(1. - 1.*Power(En,2)) + Power(L,2)/Power(x,2));
}

double AltCarterConstant(double a, double p, double e, double x){

    double En, L, Q;

    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);
    Q = CarterConstant(a,p,e,x);

    return Q + Power(L - a*En, 2);
}

double RadialRoot3(double a, double p, double e, double x){

    double En, Q, r1, r2;

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    En = Energy(a,p,e,x);
    Q = CarterConstant(a,p,e,x);

    return 1/(1. - 1.*Power(En,2)) + 0.5*(-1.*r1 - 1.*r2) + 
        Sqrt((-1.*Power(a,2)*Q)/((1. - 1.*Power(En,2))*r1*r2) + 
        Power(-1./(1. - 1.*Power(En,2)) + 0.5*(r1 + r2),2));
}

double RadialRoot4(double a, double p, double e, double x){
    double En,Q, r1,r2,r3;

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    En = Energy(a,p,e,x);
    Q = CarterConstant(a,p,e,x);
    r3 = RadialRoot3(a,p,e,x);

    return (Power(a,2)*Q)/(r1*r2*(1. - 1.*Power(En,2))*r3);
}

double PolarRoot2(double a, double p, double e, double  x){
    double En, L;

    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);

    return Sqrt(Power(a,2)*(1. - 1.*Power(En,2)) + Power(L,2)/(Power(x,2)));
}



// Mino-time Fundamental Frequencies

double RadialMinoFrequency(double a, double p, double e, double x){ 

    double r1,r2,r3,r4,En, kr;

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x); 
    En = Energy(a,p,e,x); 

    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

    return (Pi*Sqrt((1 - Power(En,2))*(r1 - r3)*(r2 - r4)))/(2.*EllipticK(kr));
}

double PolarMinoFrequency(double a, double p, double e, double x){
    
    double En, zm, zp, kz;

    zm = PolarRoot1(x);
    zp = PolarRoot2(a,p,e,x);
    En = Energy(a,p,e,x);

    kz = (Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2);

    return (Pi*zp)/(2.*EllipticK(kz));
}

double TimeMinoFrequency(double a, double p, double e, double x){

    double En, L,Q, r1, r2 ,r3 ,r4, rout, rin, zm, zp, kr, kz, hr, hout, hin, hM, 
    Kkr,Ekr, Pihrkr, PihMkr,radial,polar;

    double M = 1.0;

    rout = OuterHorizon(a);
    rin = InnerHorizon(a);
    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x);

    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

    zm = PolarRoot1(x);
    zp = PolarRoot2(a,p,e,x);

    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);
    Q = CarterConstant(a,p,e,x);
    kz = (Power(a,2)*(1.0 - Power(En,2))*Power(zm,2))/Power(zp,2);

    
    hr = (r1 - r2)/(r1 - r3);
    hout = hr * (r3 - rout)/(r2 - rout);
    hin = hr * (r3 - rin)/(r2 - rin); 
    hM =  ((r1 - r2)*(r3 - M))/((r1 - r3)*(r2 - M));

    Kkr = EllipticK(kr);
    Ekr = EllipticE(kr);
    Pihrkr = EllipticPi(hr,kr);


    radial = (4 + Power(a,2))*En + En*(2*((Pihrkr*(r2 - r3))/Kkr + r3) + 
      (-(r1*r2) + r3*(r1 + r2 + r3) + (Ekr*(r1 - r3)*(r2 - r4))/Kkr + 
         (Pihrkr*(r2 - r3)*(r1 + r2 + r3 + r4))/Kkr)/2. + 
      (-(((-2*Power(a,2) + (4 - (a*L)/En)*rin)*(1 - ((r2 - r3)*EllipticPi(hin,kr))/(Kkr*(r2 - rin))))/
            (r3 - rin)) + ((-2*Power(a,2) + (4 - (a*L)/En)*rout)*
            (1 - ((r2 - r3)*EllipticPi(hout,kr))/(Kkr*(r2 - rout))))/(r3 - rout))/Sqrt(1 - Power(a,2)));

    if(fabs(x) == 1){
        polar = - Power(a,2) * En;
    }
    else {
        polar = -(Power(a,2)*En) + (En*Q*(1 - EllipticE((Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2))/
         EllipticK((Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2))))/((1 - Power(En,2))*Power(zm,2));
    }

    return radial + polar;
}

double AzimuthalMinoFrequency(double a, double p, double e, double x){

    double En, L,Q, r1, r2 ,r3 ,r4, rout, rin, zm, zp, kr, kz, hr, hout, hin, hM, 
    Kkr,Ekr, Pihrkr, PihMkr,radial,polar;

    double M = 1.0;

    rout = OuterHorizon(a);
    rin = InnerHorizon(a);
    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x);

    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

    zm = PolarRoot1(x);
    zp = PolarRoot2(a,p,e,x);

    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);
    Q = CarterConstant(a,p,e,x);
    kz = (Power(a,2)*(1.0 - Power(En,2))*Power(zm,2))/Power(zp,2);

    
    hr = (r1 - r2)/(r1 - r3);
    hout = hr * (r3 - rout)/(r2 - rout);
    hin = hr * (r3 - rin)/(r2 - rin); 
    hM =  ((r1 - r2)*(r3 - M))/((r1 - r3)*(r2 - M));

    Kkr = EllipticK(kr);
    Ekr = EllipticE(kr);
    Pihrkr = EllipticPi(hr,kr);


    radial = (a*(-(((-(a*L) + 2*En*rin)*(1 - ((r2 - r3)*EllipticPi(hin,kr))/(Kkr*(r2 - rin))))/(r3 - rin)) + 
       ((-(a*L) + 2*En*rout)*(1 - ((r2 - r3)*EllipticPi(hout,kr))/(Kkr*(r2 - rout))))/(r3 - rout)))/
   (2.*Sqrt(1 - Power(a,2)));

   polar = (L*EllipticPi(Power(zm,2),(Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2)))/
   EllipticK((Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2));

    return radial + polar;
}


// Boyer-Linquest Time Frequencies

double RadialFrequency(double a, double p, double e, double x){ 


    return RadialMinoFrequency(a,p,e,x)/TimeMinoFrequency(a,p,e,x);
}

double PolarFrequency(double a, double p, double e, double x){ 


    return PolarMinoFrequency(a,p,e,x)/TimeMinoFrequency(a,p,e,x);
}

double AzimuthalFrequency(double a, double p, double e, double x){ 


    return AzimuthalMinoFrequency(a,p,e,x)/TimeMinoFrequency(a,p,e,x);
}
