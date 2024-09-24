double unsignedfmod(double x, double y);

double RadialMinoPhaseToDarwinPhase(double a, double p, double e, double x, double qr);

double RadialDarwinPhaseToMinoPhase(double a, double p, double e, double x, double psi);

double PolarMinoPhaseToDarwinPhase(double a, double p, double e, double x, double qz);

double PolarDarwinPhaseToMinoPhase(double a, double p, double e, double x, double chi);

double AzimuthalMinoPhaseToCoordinate(double a, double p, double e, double x, double qr,double qz, double qphi);

double AzimuthalCoordinateToMinoPhase(double a, double p, double e, double x, double psi,double chi, double phi);

double RadialMinoPhaseToBoyerLindquistPhase(double a, double p, double e, double x, double qr, double qz);

double PolarMinoPhaseToBoyerLindquistPhase(double a, double p, double e, double x, double qr, double qz);

double AzimuthalMinoPhaseToBoyerLindquistPhase(double a, double p, double e, double x, double qr, double qz,double qphi);

void DarwinPhasesToBoyerLindquistPhases(double a, double p, double e, double x, double psi, double chi,double phi, double *Phi_r, double *Phi_theta, double *Phi_phi);

void DarwinPhasesToBoyerLindquistPhasesVectorized(double* a, double* p, double* e, double* x, double* psi, double* chi,double* phi, double *Phi_r, double *Phi_theta, double *Phi_phi,int length);


void BoyerLindquistPhasesToDarwinPhases(double a, double p, double e, double x, double Phi_r, double Phi_theta,double Phi_phi, double *psi, double *chi, double *phi);

void BoyerLindquistPhasesToDarwinPhasesVectorized(double* a, double* p, double* e, double* x, double *Phi_r, double *Phi_theta, double *Phi_phi, double *psi, double *chi, double *phi, int length);


