//! Simple particle class including also previous coordinates

struct Particle{
  //!  x, y, z=s-s0
  double x, y, z;
  //! dp/p0, x', y'
  double dp, xs, ys;
  //! old x
  double x1, xs1, x2, xs2;
  //! old y
  double y1, ys1, y2, ys2;
  //! betatron wave lengths
  double lambda_h, lambda_tmp_h, lambda_v, lambda_tmp_v;
};


//! Synchronous particle

struct SynParticle{
  double Z; //!< charge number
  double A; //!< mass number
  double beta0;
  double gamma0;
  double eta0;
};   


//! Particle container

class Pic{
  vector<Particle> pics;
  double charge;  // macro particle charge/qe
  SynParticle *SP;  // synchronous particle
  FILE *out;  // output file

 public:
  //! left and right boundaries of the slice
  double z1, z2;

  //  Pic(){}  // useful?; SP
  Pic(SynParticle* ptr, double q, int n, string filename): pics(0), charge(q){
    pics.reserve(n);
    SP = ptr;
    out = fopen(filename.c_str(), "w");
  }
  Pic(SynParticle* ptr, double q, int n): pics(0), charge(q){
    pics.reserve(n);
    SP = ptr;
    out = NULL;
  }
  ~Pic(){ if(out != NULL) fclose(out); }

  double get_charge(){ return charge; }
  int get_size(){ return pics.size(); }
  vector<Particle>* get_particles(){ return &pics; }
  void clear_particles(){ pics.clear(); }


  // longitudinal distributions

  //! coasting beam with parabolic momentum distribution
  void parabolic_dc(double dummyd, double length, double dp0, long Np, long *d,
		    int dummyi=0);
  //! Parabolic bunch
  void parabolic(double zlm, double z0, double dp0, long Np, long *d, int dummyi=0);
  //! Coasting, Gaussian; SP
  void coast_gauss(double dummyd, double length, double dp0, long Np, long *d,
		   int dummyi=0);
  //! coasting beam
  void coasting_beam(double length, long Np, long *d);
  //! Gaussian momentum spread
  void gaussz(double dp0, long *d);
  //! Gaussian bunch
  void bunch_gauss(double zlm, double circum, double dp0, long Np, long *d, int dummyi=0);
  //! Constant bunch dist.
  void bunch_const(double zlm, double circum, double dp0, long Np, long *d, int linrf);
  //! barrier air bag
  void barrier_air_bag(double zlm, double dummyd, double dp0, long Np, long *d,
		       int dummyi=0);
  //! bunch air bag
  void bunch_air_bag(double zlm, double circum, double dp0, long Np, long *d, int dummyi=0);


  // Transverse distributions

  //! tranverse waterbag distribution:
  void waterbag_xy(double emittance_x, double emittance_y, double alpha_x, double alpha_y,
		   double beta_x, double beta_y, double D0, double Ds0, double centering,
		   double inj_angle, long *d);
  //! KV distribution
  void KV_xy(double emittance_x, double emittance_y, double alpha_x, double alpha_y,
	     double beta_x, double beta_y, double D0, double Ds0, double centering,
	     double inj_angle, long *d);
  //! Semi-gaussian
  void SG(double emittance_x, double emittance_y, double alpha_x, double alpha_y,
	  double beta_x, double beta_y, double D0, double Ds0, double centering,
	  double inj_angle, long *d);
  //! Gaussian
  void Gauss_xy(double emittance_x, double emittance_y, double alpha_x, double alpha_y,
		double beta_x, double beta_y, double D0, double Ds0, double centering,
		double inj_angle, long *d);

  //! basic PIC output:

  void print(int subset);

  // rms momenta:

  double rms_emittance_x();
  double rms_emittance_y();
  double x_rms();
  double x_max();
  double y_rms();
  double y_max();
  double z_mean();
  double z_min();
  double z_max();
  double z2_mean();
  double rms_z_width();
  double pz_mean();
  double pz2_mean();
  double rms_momentum_spread();
  double x2y2();
  double xy();
  double xzn(double n, double zlm);
  double offset_x();
  double offset_y();
  double entropy(Grid2D& xsys);
 
  // single particle phase advance:

  //void update_wavelength_h(double ds, double offset_x);
  //void update_wavelength_v(double ds);
  //double get_wavelength_h(int j);
  //double get_wavelength_v(int j);
  //double rms_wavelength_h();
  //double rms_wavelength_v();
 
  void init_old_coord();  // SP
  void store_old_coordinates();
  double get_phaseadvance_h(int j);
  double get_phaseadvance_v(int j);
  double rms_phaseadvance_h();
  double rms_phaseadvance_v();
 
  //! particle transport using sector map
  void transport(SectorMap* M, double K_xix, double K_xiy, double boundary);
  long localLoss_x(double loBound, double upBound);
  //! space charge kick
  void kick(Grid2D& Ex, Grid2D& Ey, double ds);
  void kick(Grid3D& Ex, Grid3D& Ey, double ds);
  //! constant kick
  void kick(double fx, double fy);
  //! cavity kick
  void cavity_kick(double voltage0, int harmonic,double R);
  void cavity_kick_linear(double voltage0, int harmonic,double R);
  //! barrier bucket kick (reflection at zm1 and zm2)
  void barrier_kick(double zm1, double zm2);
  //! thin lens kick
  void kick(ThinLens& M, double ds);
  //! transverse impedance dipole kick
  void impedance_kick(Grid1D& kick, double circum, double ds);
  //! non/linear transverse space charge kick
  void linear_SC_kick(double dQxm, double dQym, double tunex, double tuney, Grid1D& ldy,
		      double ldy0, Grid1D& dipole_current_x, Grid1D& dipole_current_y,
		      double circum, double ds);
  void nonlinear_SC_kick(double xrms, double yrms, double dQxm, double dQym, double tunex,
			 double tuney, Grid1D& ldy, double ldy0, double circum, double ds);
  void dipole_kick_simple(double dQxm, double dQym, double tunex, double tuney, Grid1D& ldy,
			  double ldy0, Grid1D& dipole_current_x, Grid1D& dipole_current_y,
			  double circum, double ds);
  //! btf noise kick
  double dipole_mod_kick(double t, double ds, double circum, double theta, double freq0,
			 double freq1, double tend, int n, long* d);
  //! btf kick
  double dipole_mod_kick(double t, double ds, double circum, double theta, double freq,
			 int n);
  //! signal from pickup 
  double pickup_signal(Grid1D& dipole_current, double circum, double t);
 
  //! add friction and diffusion (simple Fokker-Planck)
  void langevin(double beta_fxy, double beta_fz, double Dxy, double Dz, double ds,
		double betx, double bety, long* d);
 
  // grid to particle/particle to grid interpolations:

  void gatherZ(double pic_charge,Grid1D& target);
  void gatherX(double pic_charge,Grid1D& target);
  void gatherXs(double pic_charge,Grid1D& target);
  void gatherY(double pic_charge,Grid1D& target);
  void gatherXY(double pic_charge,Grid2D& target);
  void gatherXYZ(double pic_charge,Grid3D& target);
  void gatherXXs(double pic_charge,Grid2D& target);
  void gatherZX(double pic_charge,Grid2D& target);
  void gatherYYs(double pic_charge,Grid2D& target);
  void gatherXsYs(double pic_charge,Grid2D& target);

  // exchange of particles

  vector<Particle> get_particles_left(double length);
  vector<Particle> get_particles_right(double length);
  void add_particles(vector<Particle>& part);
  void remove_lost_particles(double length);
  void periodic_bc(double length);
};
