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
        printf("Usage: ./KerrCoM a p e x psi chi phi\n");
    }

    clock_t start, end;
    double a, p, e, x, r1, r2, r3, r4, zm, cpu_time_used, 
    qr,qz,qphi,Freq_r,Freq_z,r,z,theta,t,Freq_t,Freq_phi, 
    psi,chi,phi,test,Phi_r,Phi_theta,Phi_phi;

    a = atof(argv[1]);
    p = atof(argv[2]);
    e = atof(argv[3]);
    x = atof(argv[4]);
    psi = atof(argv[5]);
    chi = atof(argv[6]);
    phi = atof(argv[7]);

    printf("___Parameters Entered___\n");
    printf("Spin: %f\n", a);
    printf("Semi-Lattice Rectum: %f\n", p);
    printf("Eccentricity: %f\n", e);
    printf("Cos(Ang_Inc) %f\n", x);
    printf("Radial Phase: %f\n", psi);
    printf("Polar Phase %f\n", chi);
    printf("Azimuthal Phase %f\n", phi);



    printf("\n___Angle Conversions___\n");
    start = clock();
    DarwinPhasesToBoyerLindquistPhases(a,p,e,x,psi,chi,phi,&Phi_r,&Phi_theta,&Phi_phi);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC *1000;
    printf("Radial Boyer-Lindquist Phase: %f\n", Phi_r);
    printf("Radial Boyer-Lindquist Phase: %f\n", Phi_theta);
    printf("Radial Boyer-Lindquist Phase: %f\n", Phi_phi);
    printf("CPU time taken: %f ms\n", cpu_time_used);
    
    return 0;
}
