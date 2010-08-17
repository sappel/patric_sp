//! Induced xprime kick per turn
void InducedWakeKick(Grid1D& kick,Grid1D& dipole_current,Grid1D& dipole_current_x,
                 double tune, double omega0, double nres, double Rs, double Qs, double b, double leit, double beta0, double E0, double charge); 
void InducedKick(Grid1D& kick,Grid1D& dipole_current,Grid1D& dipole_current_s,
                 double tune, double omega0, double nres, double Rs, double Qs, double beta0, double E0, double charge, double t); 
void InducedKick(Grid1D& kick,Grid1D& dipole_current,
                 double tune, double omega0, double nres, double Rs, double Qs, double beta0, double E0, double charge); 
void InducedKick(Grid1D& kick,Grid1D& dipole_current,
                 double Zimage, double beta0, double E0, double charge); 
//! Induced xprime kick per turn: coasting beam approximation
double InducedKick(double offset, double ds, komplex dqc,double beta0, double Qx, double circum);


