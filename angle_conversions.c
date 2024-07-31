#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <math.h>
#include<float.h>
#include "common_definitions.h"
#include "constants_of_motion.h"
#include "oscillating_terms.h"

double unsignedfmod(double x, double y){
    if(x>= 0){
        return fmod(x,y);
    }
    else{ 
        return 2*Pi + fmod(x,y);
    }
}

double RadialMinoTimePhaseToQuasiKeplerian(double a, double p, double e, double x, double qr){ 

    double bulk,r, remainder, psi;
    
    bulk = qr - unsignedfmod(qr , 2*Pi);

    r = RadialMinoTimeSol(a,p,e,x,qr);

    psi = acos((p/r -1)/e);

    if(unsignedfmod(qr , 2*Pi) - Pi <= 0){
        remainder = psi;
    }
    else{
        remainder = 2*Pi - psi;
    }

    return bulk + remainder;
    }

double RadialQuasiKeplerianPhaseToMinoTime(double a, double p, double e, double x, double psi){ 

    double bulk,r, remainder, qr,r1,r2,r3,r4,kr,kpsi;

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x); 

    kr = ((r1 - r2)*(r3 - r4))/((r1 - r3)*(r2 - r4));

    r = RadialQuasiKeplerianSol(a,p,e,x,psi);

    kpsi = Sqrt( ((r - r2)*(r1 - r3))/((r - r3)*(r1 - r2)));

    qr = (Pi/EllipticK(kr)) * EllipticF(asin(kpsi),kr);

    bulk = psi - unsignedfmod(psi , 2*Pi);

    if(unsignedfmod(psi , 2*Pi) - Pi <= 0){
        remainder = qr;
    }
    else{
        remainder = 2*Pi - qr;
    }

    return bulk + remainder;
    }


double PolarMinoTimePhaseToQuasiKeplerian(double a, double p, double e, double x, double qz){ 

    double bulk,z,zm, remainder, chi;
    
    bulk = qz - unsignedfmod(qz , 2*Pi);
    zm = PolarRoot1(x);

    z = PolarMinoTimeSol(a,p,e,x,qz); 

    chi= acos(z/zm); 

    if(unsignedfmod(qz , 2*Pi) - Pi <= 0){
        remainder = chi;
    }
    else{
        remainder = 2*Pi - chi;
    }

    return bulk + remainder;
    }

double PolarQuasiKeplerianPhaseToMinoTime(double a, double p, double e, double x, double chi){ 

    double bulk,z, remainder, qz,zm,zp, En, kz;

    zm = PolarRoot1(x); 
    zp = PolarRoot2(a,p,e,x);
    En = Energy(a,p,e,x);

    kz = (Power(a,2)*(1 - Power(En,2))*Power(zm,2))/Power(zp,2);

    z = PolarQuasiKeplerianSol(a,p,e,x,chi);

    qz = - ( (Pi * EllipticF(asin(cos(chi)),kz))/(2 * EllipticK(kz)) - Pi/2);

    bulk = chi - unsignedfmod(chi , 2*Pi);

    if(unsignedfmod(chi , 2*Pi) - Pi <= 0){
        remainder = qz;
    }
    else{
        remainder = 2*Pi - qz;
    }

    return bulk + remainder;
    }

double AzimuthalMinoTimePhaseToQuasiKeplerian(double a, double p, double e, double x, double qr,double qz, double qphi){

     // Return phi coordinate
    return qphi + DeltaPhi(a,p,e,x,qr,qz);
   }

double AzimuthalQuasiKeplerianPhaseToMinoTime(double a, double p, double e, double x, double psi,double chi, double phi){

    double qr, qz,qphi;

    qr = RadialQuasiKeplerianPhaseToMinoTime(a,p,e,x,psi);
    qz = PolarQuasiKeplerianPhaseToMinoTime(a,p,e,x,chi);

    qphi = phi - DeltaPhi(a,p,e,x,qr,qz);

    return qphi;
   }

