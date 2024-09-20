# Declare external function from mylib
cdef extern from "../include/angle_conversions.h":
     void DarwinPhasesToBoyerLindquistPhases(double a, double p, double e, double x, double psi, double chi,double phi, 
        double *Phi_r, double *Phi_theta, double *Phi_phi);


# Python Wrapper
def pyDarwinPhasesToBoyerLindquistPhases(double a, double p, double e, double x, double psi, double chi,double phi):
    cdef double Phi_r
    cdef double Phi_theta
    cdef double Phi_phi

    DarwinPhasesToBoyerLindquistPhases(a,p,e,x,psi,chi,phi, &Phi_r, &Phi_theta, &Phi_phi)

    return (Phi_r,Phi_theta,Phi_phi)