#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <math.h>
#include<float.h>
#include "constants_of_motion.h"

// Definitions for using CForm'ed output from Mathematica
#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Pi				M_PI


// Radial Functions

double RadialMinoTimeSol(double a, double p, double e, double x, double qr){ 

    double r1,r2,r3,r4,kr, JacobiSN, JacobiCN, JacobiDN, u;

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x); 

    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));
    u = (qr*EllipticK(kr))/Pi;
    gsl_sf_elljac_e(u, kr, &JacobiSN, &JacobiCN, &JacobiDN);

    return (-(r2*(r1 - r3)) + (r1 - r2)*r3*Power(JacobiSN,2))/
    (-r1 + r3 + (r1 - r2)*Power(JacobiSN,2));
}

double RadialQuasiKeplerianSol(double a, double p, double e, double x, double psi){

    return p/(1 + e * cos(psi));
}

// Polar Terms

double PolarMinoTimeSol(double a, double p, double e, double x, double qz){
    
    double En, zm, zp, kz, u, JacobiSN, JacobiCN, JacobiDN;

    zm = PolarRoot1(x);
    zp = PolarRoot2(a,p,e,x);
    En = Energy(a,p,e,x);

    kz = (Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2);
    u = (2*(qz+Pi/2)*EllipticK(kz))/Pi;

    gsl_sf_elljac_e(u, kz, &JacobiSN, &JacobiCN, &JacobiDN);


    return zm*JacobiSN;
}

double PolarQuasiKeplerianSol(double a, double p, double e, double x, double chi){

    double zm = PolarRoot1(x);

    return zm*cos(chi);
}


// Time Coordiante Quantities
double DeltaT(double a, double p, double e, double x, double qr, double qz){

    double En, L, r1, r2 ,r3 ,r4, rp, rm, zm, zp, kr, kz, hr, hp, hm, 
    UpsilonZ, UpsilonR, K_kr, E_kr, Pi_hr_kr, Pi_hM_kr, Gamma, 
    snR, snZ, cnR, dnR, cnZ, dnZ, psiR, psiZ, tr_qr, tz_qz, ur, uz;

    double M = 1.0;
    rp = OuterHorizon(a);
    rm = InnerHorizon(a);
    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x);
    zm = PolarRoot1(x);
    zp = PolarRoot2(a,p,e,x);
    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);
    UpsilonR = RadialMinoFrequency(a,p,e,x);
    UpsilonZ = PolarMinoFrequency(a,p,e,x);
    Gamma = TimeMinoFrequency(a,p,e,x);


    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

    kz = (Power(a,2)*(1.0 - Power(En,2))*Power(zm,2))/Power(zp,2);

    hr = (r1 - r2)/(r1 - r3);
    hp = hr * (r3 - rp)/(r2 - rp);
    hm = hr * (r3 - rm)/(r2 - rm);

    // Calculating Jacobi Ampltudes
    ur = EllipticK(kr)/Pi *qr;
    uz = (2*(qz + Pi/2)*EllipticK(kz))/Pi;
    
   psiR = Jacobi_am(ur,'m',kr);
   psiZ = Jacobi_am(uz,'m',kz);
    

    // Radial and Polar Oscilations
    tr_qr = -((En*((-4*(r2 - r3)*(-(((-2*Power(a,2) + (4 - (a*L)/En)*rm)*
                   ((qr*EllipticPi(hm,kr))/Pi - 
                     EllipticPiIncomp(hm,psiR,kr)))/
                 ((r2 - rm)*(r3 - rm))) + 
              ((-2*Power(a,2) + (4 - (a*L)/En)*rp)*
                 ((qr*EllipticPi(hp,kr))/Pi - 
                   EllipticPiIncomp(hp,psiR,kr)))/((r2 - rp)*(r3 - rp)))
            )/(-rm + rp) + 
         4*(r2 - r3)*((qr*EllipticPi(hr,kr))/Pi - 
            EllipticPiIncomp(hr,psiR,kr)) + 
         (r2 - r3)*(r1 + r2 + r3 + r4)*
          ((qr*EllipticPi(hr,kr))/Pi - EllipticPiIncomp(hr,psiR,kr)) + 
         (r1 - r3)*(r2 - r4)*
          ((qr*EllipticE(kr))/Pi - EllipticEIncomp(psiR,kr) + 
            (hr*Cos(psiR)*Sin(psiR)*Sqrt(1 - kr*Power(Sin(psiR),2)))/
             (1 - hr*Power(Sin(psiR),2)))))/
     Sqrt((1 - Power(En,2))*(r1 - r3)*(r2 - r4)));

    tz_qz = (En*zp*((2*(qz+ Pi/2)*EllipticE(kz))/Pi - EllipticEIncomp(psiZ,kz)))/
   (1 - Power(En,2));

    
    // Geodesic Solution
    return tr_qr + tz_qz;
}

// Azimuthal Functions

double DeltaPhi(double a, double p, double e, double x, double qr, double qz){

    double En, L, r1, r2 ,r3 ,r4, rp, rm, zm, zp, kr, kz, hr, hp, hm, 
    UpsilonZ, UpsilonR, K_kr, E_kr, Pi_hr_kr, Pi_hM_kr, UpsilonPhi, 
    snR, snZ, cnR, dnR, cnZ, dnZ, psiR, psiZ, phir_qr, phiz_qz, ur, uz;

    double M = 1.0;
    rp = OuterHorizon(a);
    rm = InnerHorizon(a);
    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x);
    zm = PolarRoot1(x);
    zp = PolarRoot2(a,p,e,x);
    En = Energy(a,p,e,x);
    L = AngularMomentum(a,p,e,x);
    UpsilonR = RadialMinoFrequency(a,p,e,x);
    UpsilonZ = PolarMinoFrequency(a,p,e,x);
    UpsilonPhi = AzimuthalMinoFrequency(a,p,e,x);


    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

    kz = (Power(a,2)*(1.0 - Power(En,2))*Power(zm,2))/Power(zp,2);

    hr = (r1 - r2)/(r1 - r3);
    hp = hr * (r3 - rp)/(r2 - rp);
    hm = hr * (r3 - rm)/(r2 - rm);

    // Calculating Jacobi Ampltudes
    ur = EllipticK(kr)/Pi *qr;
    uz = (2*(qz + Pi/2)*EllipticK(kz))/Pi;
    
   psiR = Jacobi_am(ur,'m',kr);
   psiZ = Jacobi_am(uz,'m',kz);
    

    // Radial and Polar Oscilations
    phir_qr = (2*a*En*(-(((r2 - r3)*(-((a*L)/En) + 2*rm)*
            ((qr*EllipticPi(hm,kr))/Pi - EllipticPiIncomp(hm,psiR,kr)))/
          ((r2 - rm)*(r3 - rm))) + 
       ((r2 - r3)*(-((a*L)/En) + 2*rp)*
          ((qr*EllipticPi(hp,kr))/Pi - EllipticPiIncomp(hp,psiR,kr)))/
        ((r2 - rp)*(r3 - rp))))/
   (Sqrt((1 - Power(En,2))*(r1 - r3)*(r2 - r4))*(-rm + rp));

    phiz_qz = -((L*((2*(qz +Pi/2)*EllipticPi(Power(zm,2),kz))/Pi - 
         EllipticPiIncomp(Power(zm,2),psiZ,kz)))/zp);

    
    // Geodesic Solution
    return phir_qr + phiz_qz;
}