#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common_definitions.h"
#include "constants_of_motion.h"
#include "oscillating_terms.h"
#include "angle_conversions.h"
#include <time.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>



int main(int argc, char **argv){

    if(argc < 7){
        printf("Error: Not enough arguments.\n");
        printf("Usage: ./KerrCoM a p e x q_r q_z q_phi\n");
    }

    clock_t start, end;
    double a, p, e, x, r1, r2, r3, r4, zm, cpu_time_used, 
    qr,qz,qphi,Freq_r,Freq_z,r,z,theta,t,Freq_t,Freq_phi, 
    psi,chi,phi,test;

    a = atof(argv[1]);
    p = atof(argv[2]);
    e = atof(argv[3]);
    x = atof(argv[4]);
    qr = atof(argv[5]);
    qz = atof(argv[6]);
    qphi = atof(argv[7]);

    printf("___Parameters Entered___\n");
    printf("Spin: %f\n", a);
    printf("Semi-Lattice Rectum: %f\n", p);
    printf("Eccentricity: %f\n", e);
    printf("Cos(Ang_Inc) %f\n", x);
    printf("Radial Phase: %f\n", qr);
    printf("Polar Phase %f\n", qz);
    printf("Azimuthal Phase %f\n", qphi);

    r1 = RadialRoot1(p,e);
    r2 = RadialRoot2(p,e);
    r3 = RadialRoot3(a,p,e,x);
    r4 = RadialRoot4(a,p,e,x);
    zm = PolarRoot1(x);
    Freq_r = RadialMinoFrequency(a,p,e,x);
    Freq_z = PolarMinoFrequency(a,p,e,x);
    Freq_t = TimeMinoFrequency(a,p,e,x);
    Freq_phi = AzimuthalMinoFrequency(a,p,e,x);

    r = RadialMinoTimeSol(a,p,e,x,qr);
    z = PolarMinoTimeSol(a,p,e,x,qz);
    theta = acos(z);

    /*
    printf("\n___Frequencies___\n");
    printf("Radial: %f\n", Freq_r);
    printf("Polar:  %f\n", Freq_z);
    printf("Time:  %f\n", Freq_t);
    printf("Azimuthal:  %f\n", Freq_phi);
    
    printf("\n___Analytic Solutions___\n");
    printf("t = %f\n", t);
    printf("r =  %f\n", r);
    printf("theta =  %f\n", theta);
    printf("Phi = %f\n", phi);
    */
    printf("\n___Angle Conversions___\n");

    psi = RadialMinoTimePhaseToQuasiKeplerian(a,p, e,x,qr);
    printf("psi = %f\n", psi);
    chi = PolarMinoTimePhaseToQuasiKeplerian(a,p, e,x,qz);
    printf("chi = %f\n", chi); 
    test = RadialQuasiKeplerianPhaseToMinoTime(a,p,e,x,psi);
    printf("qr = %f\n", test);
    printf("Differnece = %f\n", test - qr);

    test = PolarQuasiKeplerianPhaseToMinoTime(a,p,e,x,chi);
    printf("qz = %f\n", test);
    printf("Differnece = %f\n", test - qz);

    phi = AzimuthalMinoTimePhaseToQuasiKeplerian(a,p,e,x,qr,qz,qphi);
    printf("phi = %f\n", phi);
    test = AzimuthalQuasiKeplerianPhaseToMinoTime(a,p,e,x,psi,chi,phi);
    printf("qphi = %f\n", test);
     printf("Differnece = %f\n", test - qphi);
    return 0;
}
