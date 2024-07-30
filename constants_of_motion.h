double Jacobi_am(double u, char arg,  double x);

//Horizon Locations

double InnerHorizon(double a);

double OuterHorizon(double a);

//Primary Roots

double RadialRoot1(double p, double e);

double RadialRoot2(double p, double e);

double PolarRoot1(double x);

double Delta(double a, double r);

// Functions needed for E,L,Q calculations

double D(double a, double r, double zm);

double F(double a, double r, double zm);

double G(double a, double r);

double H(double a, double r, double zm);

// E, L, Q Calcuaitons

double Energy(double a, double p, double e, double x);

double AngularMomentum(double a, double p, double e, double x);

double CarterConstant(double a, double p, double e, double x);

double AltCarterConstant(double a, double p, double e, double x);

double RadialRoot3(double a, double p, double e, double x);

double RadialRoot4(double a, double p, double e, double x);

double PolarRoot2(double a, double p, double e, double x);

// Elliptic Functions

double EllipticK(double k);

double EllipticF(double phi, double k);

double EllipticE(double k);

double EllipticEIncomp(double phi, double k);

double EllipticPi(double n, double k);

double EllipticPiIncomp(double n, double phi, double k);

// Mino Time Frequencies

double TimeMinoFrequency(double a, double p, double e, double x);

double RadialMinoFrequency(double a, double p, double e, double x);

double PolarMinoFrequency(double a, double p, double e, double x);

double AzimuthalMinoFrequency(double a, double p, double e, double x);

// Boyer-Linquist Time Frequencies

double RadialFrequency(double a, double p, double e, double x);

double PolarFrequency(double a, double p, double e, double x); 

double AzimuthalFrequency(double a, double p, double e, double x);