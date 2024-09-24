import numpy as np
cimport numpy as np
# from libcpp.string cimport string
# from libcpp cimport bool

assert sizeof(int) == sizeof(np.int32_t)


# Declare external function from mylib
cdef extern from "../include/angle_conversions.h":
     void DarwinPhasesToBoyerLindquistPhases(double a, double p, double e, double x, double psi, double chi,double phi, 
        double *Phi_r, double *Phi_theta, double *Phi_phi);

     void BoyerLindquistPhasesToDarwinPhases(double a, double p, double e, double x, double Phi_r, double Phi_theta,double Phi_phi, 
        double *psi, double *chi, double *phi);


# Python Wrapper
def pyDarwinPhasesToBoyerLindquistPhases(double a, 
                                         double p, 
                                         double e, 
                                         double x, 
                                         double psi, 
                                         double chi,
                                         double phi):

    cdef double Phi_r
    cdef double Phi_theta 
    cdef double Phi_phi 

    DarwinPhasesToBoyerLindquistPhases(a,p,e,x,psi,chi,phi, &Phi_r, &Phi_theta, &Phi_phi)

    return (Phi_r,Phi_theta,Phi_phi)

def pyBoyerLindquistPhasesToDarwinPhases(double a, double p, double e, double x,double Phi_r, double Phi_theta,double Phi_phi):
    cdef double psi
    cdef double chi
    cdef double phi

    BoyerLindquistPhasesToDarwinPhases(a,p,e,x,Phi_r,Phi_theta,Phi_phi, &psi, &chi, &phi)

    return (psi,chi,phi)

# Declare external function from mylib
cdef extern from "../include/oscillating_terms.h":
    double RadialDarwinFrequency(double a, double p, double e, double x, double psi, double chi);
    double PolarDarwinFrequency(double a, double p, double e, double x, double psi, double chi);
    double AzimuthalDarwinFrequency(double a, double p, double e, double x, double psi, double chi);

def pyRadialDarwinFrequency(double a, double p, double e, double x,double psi, double chi):

    return RadialDarwinFrequency(a,p,e,x,psi,chi)
    
def pyPolarDarwinFrequency(double a, double p, double e, double x,double psi, double chi):

    return PolarDarwinFrequency(a,p,e,x,psi,chi)

def pyAzimuthalDarwinFrequency(double a, double p, double e, double x,double psi, double chi):

    return AzimuthalDarwinFrequency(a,p,e,x,psi,chi)

# Declare external function from mylib
cdef extern from "../include/constants_of_motion.h":
    double RadialFrequency(double a, double p, double e, double x);
    double PolarFrequency(double a, double p, double e, double x);
    double AzimuthalFrequency(double a, double p, double e, double x);

def pyRadialFrequency(double a, double p, double e, double x):

    return RadialFrequency(a,p,e,x)

def pyPolarFrequency(double a, double p, double e, double x):

    return PolarFrequency(a,p,e,x)

def pyAzimuthalFrequency(double a, double p, double e, double x):

    return AzimuthalFrequency(a,p,e,x)