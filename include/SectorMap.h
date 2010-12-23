//! Twiss parameters for a lattice element

struct TwissP{
  double betx, bety, alpx, alpy, Dx;  //!< beta_x, beta_y, alpha_x, alpha_y, Dispersion_x
  double mux, muy, xix, xiy;  //!< phase advance and chromaticity; SP
};

//! MAD-like Sectormap objects

class SectorMap{
  string ElementName;  //!< element name
  double L;  //!< element length
  double T[36];  //!< 1th order transport matrix
  double K[6];  //!< kick vector
  TwissP twiss;  //!< Twiss parameter

 public:
  //! construct empty map
  SectorMap(){
    L = 0.0;
    for(int j=0; j<36; j++) T[j] = 0.0;
    for(int j=0; j<6; j++) K[j] = 0.0;
  }
  //! copy constructor
  SectorMap(const SectorMap& S){
    ElementName = S.ElementName;	
    L = S.L;
    for(int j = 0; j<36; j++) T[j] = S.T[j];
    for(int j = 0; j<6; j++) K[j] = S.K[j];
    twiss = S.twiss;
  }
  //! construct a simple continuous focusing element
  SectorMap(double sigx, double sigy, double R, double length, double gamma0);
  ~SectorMap(){};

  SectorMap& operator=(const SectorMap& M);
  //! multiply two maps
  SectorMap operator*(const SectorMap& M);
  //! get or set element name
  string& get_name(){ return ElementName; }
  //! get or set element length
  double& get_L(){ return L; }
  double& get_K(int j){ return K[j]; }
  //! get or set matrix elements
  double& get_T(int j, int i){ return T[j*6+i]; }
  TwissP& get_twiss(){ return twiss; }
  double get_betx(){ return twiss.betx; }
  double get_bety(){ return twiss.bety; }
  //! return pointer to sectormap
  SectorMap* get_map(){ return this; }
  //! calculate phase advances per element
  void phase_advance(double& sigx, double& sigy);
  //! transport particle coordinates through the sector
  void transport(vektor& R1, vektor& R0, double K_xix, double K_xiy);

  // the following functions have been added by SP
  double get_mux() const { return twiss.mux; }
  double get_muy() const { return twiss.muy; }
  double get_xix() const { return twiss.xix; }
  double get_xiy() const { return twiss.xiy; }
  double get_alpx() const { return twiss.alpx; }
};


//! Container class (List) for sector maps

class BeamLine{
  list<SectorMap> line;
  list<SectorMap>::iterator element;

 public:
  BeamLine(){}
  //  BeamLine(string dir);  //!< generate beam line from MADX files
  BeamLine(string dir, double &circum, double &Q_hor);  //!< generate beam line from MADX files
  BeamLine(const BeamLine& B){ line = B.line; element = B.element; }
  ~BeamLine(){}

  void init(BeamLine& B){ line = B.line; element = B.element; }
  void init(string dir, double &circum, double &Q_hor, double &Q_ver, double beta, int chroma);
  void read_madx_twiss(string fname, double &circum, double &Q_hor, double &Q_ver, double beta, int chroma);  //!< read madx twiss file
  void read_madx_sectormap(string fname);  //!< read madx sectormap file
  int get_size(){ return line.size(); }
  double get_L();
  list<SectorMap>::iterator get_element(){ return element; }
  list<SectorMap>::iterator get_first_element(){ return line.begin(); }
  list<SectorMap>::iterator get_end_element(){ return line.end(); }
  void next_element();
  void first_element(){ element = line.begin(); }
  void last_element(){ element = line.end(); --element; }
  void set_element(int j);
  void add_map(SectorMap& M){ line.push_back(M); }  //!< add a sectormap

  //SectorMap submap(int start, int end);  //!< multiply sector maps
  void phase_advance(double& sigx, double& sigy);

  void print_twiss();
  void print_T();
};


//! Thin lens elements (linear and nonlinear)

class ThinLens{
 public:	
  virtual void kick(vektor& R1, vektor& R0, double ds) = 0;
};	


class Octupole: public ThinLens{
  double strength_h;

 public:
  Octupole(double strength){ strength_h = strength; }
  void kick(vektor& R1, vektor& R0, double ds);
};


class Chrom: public ThinLens{
  double chrom_x, chrom_y; 	

 public:
  Chrom(double tunex0, double tuney0, double R);
  void kick(vektor& R1, vektor& R0, double ds);
};


class TuneShift: public ThinLens{
  double coeffx, coeffy; 	

 public:
  TuneShift(double tunex0, double tuney0, double dtunex, double dtuney,
	    double R);
  void kick(vektor& R1, vektor& R0, double ds);
};	


class AmplitudeDetuning: public ThinLens{
  double emitx, emity, coeffx, coeffy;
  SectorMap SMap;

 public:
  AmplitudeDetuning(double tunex0, double tuney0, double Qxs, double Qys,
		    double R, SectorMap& M);
  void kick(vektor& R1, vektor& R0, double ds);
};	

