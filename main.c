#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants_of_motion.h"
#include "OscillatingTerms.h"
#include <time.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>


#define Pi				M_PI


int main(int argc, char **argv){

    if(argc < 6){
        printf("Error: Not enough arguments.\n");
        printf("Usage: ./KerrCoM a p e x q_r q_z\n");
    }

    clock_t start, end;
    double a, p, e, x, r1, r2, r3, r4, zm, cpu_time_used, 
    qr,qz, Freq_r,Freq_z, r,z,theta,phi,t,Freq_t,Freq_phi;

    a = atof(argv[1]);
    p = atof(argv[2]);
    e = atof(argv[3]);
    x = atof(argv[4]);
    qr = atof(argv[5]);
    qz = atof(argv[6]);

    printf("___Parameters Entered___\n");
    printf("Spin: %f\n", a);
    printf("Semi-Lattice Rectum: %f\n", p);
    printf("Eccentricity: %f\n", e);
    printf("Cos(Ang_Inc) %f\n", x);
    printf("Radial Phase: %f\n", qr);
    printf("Polar Phase %f\n", qz);

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
    t = DeltaT(a,p,e,x,qr,qz);
    phi = DeltaPhi(a,p,e,x,qr,qz);


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
    

    
    return 0;
}
