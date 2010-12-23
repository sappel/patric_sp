#include "mpi.h"
#include "obflib.h"
#include "SectorMap.h"
#include "Pic.h"
#include "TImpedance.h"
#include "Bump.h"

#include <stdlib.h>

const double rp = qe*qe/(4*PI*eps0*mp*clight*clight);

//-------- Input variables:---------------------

/*!
  \defgroup Particles Global variables: particles/fields
*/
//@{
int NPIC;  //!< Total number of particles
int NX;  //!< Grid points in x, used by 2D poisson solver
int NY;  //!< Grid points in y, used by 2D poisson solver
int NZ;  //!< Grid points in z, used for dipole current and impedances
int NZ_bunch = NZ/2;  //!< number of slices along the bunch for 2.5D space charge solver (int multiple of numprocs )
//@}

/*!
  \defgroup Lattice Global variables: lattice/ring
*/
//@{
double piperadius;  //!< radius of the pipe [m]
double coll_halfgap;  //!< half gap width of collimator [m]
double image_x;  //!< horizontal boundary for Green's function; 0 means open
double image_y;  //!< vertical boundary for Green's function; 0 means open
double circum;  //!< ring circumference [m]
double gamma_t;  //!< gamma transition
double CF_advance_h, CF_advance_v;  //!< const. focusing phase advance per cell [rad] (madx_input_file = 0)
double CF_R = 0.0;  //!< radius for dispersion > 0  (madx_input_file = 0)
double CF_length;  //!< length of the constant focusing cell (madx_input_file = 0)
int NCF = 32;  //!< number of constant focusing sub-cells (madx_input_file = 0)
double koct = 8.0;  //!< octupole strength (octupole_kick = 1)
double dQxm = -0.15;  //!< dc beam maximum linear sc tune shift (space_charge = 2, 3)
double dQym = -0.15;  //!< dc beam maximum linear sc tune shift (space_charge = 2, 3)
double dqx_detune = 0.04, dqy_detune = 0.04;  //!< rms amplitude detuning (ampdetun_kick = 1)
//@}

/*!
  \defgroup Simulation Global variables: Output/Simulation
*/
//@{
char ausgabe[50];  //!< output directory
int pic_subset;  //!< numper of particles for output
int cells;  //!< length of the simulation run: number of cells
int print_cell = 0;  //!< output of particles every nth cell
double lossTol = 0.;  //!< tolerable relative losses
//@}

/*!
  \defgroup Beam Global variables: Beam
*/
//@{
double e_kin;  //!< kinetic energy in MeV/u
double Z;  //!< charge state
double A;  //!< mass number
double current;  //!< dc beam current in Amps.
int init_pic_xy;  //!< type of transverse beam distribution
int init_pic_z;  //!< type of longitudinal beam distribution
double momentum_spread;  //!< rms momentum spread
double rms_emittance_x0, rms_emittance_y0;  //!< rms emittances
double mismatch_x, mismatch_y;  //!< transverse mismatch
double x_septum = 0.0;  //!< distance of septum from nominal orbit [m]
double offcenter = 0.0;  //!< offset with respect to ideal injection [m]
double inj_angle = 0.0;  //!< injection angle
unsigned max_inj = 0;  //!< maximal number of injections
double inj_phase = 0.;  //!< Phase of injected beamletts in xx'-space
unsigned sept_return = 0;  //!< Number of revoultions till first return to septum (close to 1/Q_f)
double bunchfactor = 1.0;  //!< bunching factor
//@}

/*!
  \defgroup Impedance Global parameters: transverse impedance
*/
//@{
double dqci = 0.0;  //!< (external) imaginary coherent tune shift (sliced = 0)
double dqcr =-0.15;  //!< (external) real coherent tune shift (sliced = 0)
double Rs = 1.0e6;  //!< broadband oscillator peak impedance (sliced = 1)
double nres = 10.0;  //!< harmonic number of the oscillator resonant frequency (sliced = 1)
double Qs = 1.0;  //!< quality factor of the oscillator (sliced = 1)
double Zimage = 0.0;  //!< imaginary part of the image current impedance (sliced = 1)
double leit = 1.0e6;  //!< pipe wall conductivity (sliced = 1)
double dwall = 0.0003;  //!< pipe wall thickness  (sliced = 1)
//@}

/*!
  \defgroup Other Global variables: other simulation controls
*/
//@{
int madx_input_file = 0;  //!< madx input file (yes:1, no:0)
int space_charge = 3;  //!< space charge yes(1)/no(0)/linear(2)/nonlinear(3)
int imp_kick = 1;  //!< impedance on (1)
int sliced = 0;  //!< 2D (0) or 3D (1)
int cavity = 0;  //!< rf cavity(1)/no(0)/barrier(2)
int octupole_kick = 1;  //!< octupole kick off (0), off (1)
int ampdetun_kick = 0;  //!< amplitude detuning kick 0/1
int chroma = 0;  //!< chromaticity correction 0/1 (madx_input_file = 0)
int bc_end = 1;  //!< boundary condition in z: bunch(0)/periodic(1)
int footprint = 0;  //!< scheme: map(0)/wavelength(1)
//@}

/*!
  \defgroup Diagnostics variables: Pickups, BTF, etc.
*/
//@{	
int btf = 0;  //!< btf dipole noise excitation off (0), on (1)
int btf_harmonic = 0;  //!< btf mode number
//@}	


//-------------------end input ---------------------------

// loading input from cfg file (e.g. patric.cfg):

void input_from_file(string filename, int id){
  char dummy_string[80];
  FILE *cfg_file_ptr;
  long i, condition;
  string cfg_filename = filename + ".cfg";
  cfg_file_ptr = fopen (cfg_filename.c_str(), "r");
  if(cfg_file_ptr == NULL){
    printf("\n\n    Configuration-file does not exist.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &NPIC);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &NX);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &NY);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &NZ);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &NZ_bunch) ;
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &cells);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &lossTol);  // introduced by SP
  lossTol = 1 - lossTol;
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &e_kin);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &Z);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &A);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &current);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &piperadius);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &coll_halfgap);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &image_x);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &image_y);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &circum);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &gamma_t);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &CF_advance_h);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &CF_advance_v);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &CF_R);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &CF_length);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &NCF);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &koct);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &dQxm);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &dQym);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &dqx_detune);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &dqy_detune);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &pic_subset);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &init_pic_xy);
  fscanf(cfg_file_ptr, "%s", dummy_string);  	
  fscanf(cfg_file_ptr, "%d", &init_pic_z);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &momentum_spread);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &rms_emittance_x0);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &rms_emittance_y0);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &mismatch_x);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &mismatch_y);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &x_septum);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &offcenter);  // meaning changed by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &inj_angle);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%u", &max_inj);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &inj_phase);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%u", &sept_return);  // introduced by SP
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &bunchfactor);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &dqci);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &dqcr);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &Rs);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &nres);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &Qs);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &leit);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%lf", &Zimage);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &madx_input_file);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &space_charge);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &imp_kick);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &sliced);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &cavity);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &octupole_kick);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &ampdetun_kick);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &chroma);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &bc_end);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &print_cell);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &footprint);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &btf);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%d", &btf_harmonic);
  fscanf(cfg_file_ptr, "%s", dummy_string);
  fscanf(cfg_file_ptr, "%s", ausgabe);
  if(id == 0)
    cout << "Ausgabedatei:" << ausgabe << endl;
  fclose (cfg_file_ptr);
}

// write an output info file idl.dat for IDL plotting

void print_IDL(string data_dir, int numprocs, double cell_length, int Nelements){
  string idl_data = data_dir + "idl.dat";
  FILE *out = fopen(idl_data.c_str(), "w");

  fprintf(out, "%d\n", numprocs);
  fprintf(out, "%g\n", e_kin);
  fprintf(out, "%g\n", Z);
  fprintf(out, "%g\n", A);
  fprintf(out, "%g\n", current);
  fprintf(out, "%g\n", circum);
  fprintf(out, "%d\n", Nelements);
  fprintf(out, "%g\n", cell_length);
  fprintf(out, "%d\n", pic_subset);
  fprintf(out, "%g\n", piperadius);
  fprintf(out, "%d\n", NX);
  fprintf(out, "%d\n", NY);
  fprintf(out, "%d\n", NZ);
  fprintf(out, "%d\n", print_cell);

  fflush(out);
  fclose(out);
}


int gdbflag;
void waitforgdb(int myid){
  //myid=1;  // catch all processes to let the 0th crash before the others; execution won't finish
  if(myid == 0){
    printf("PID %d ready for attaching.\n", getpid());
    fflush(stdout);
    sync();
    gdbflag = 1;  // 1: activated, else: deactivated
    while(1 == gdbflag)
      sleep(1);
  }
}


// main part of the Simulation:

main(int argc, char* argv[]){
  time_t time1 = time(0), time2;

  //-------MPI initialzation-------------

  int numprocs, myid, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  MPI_Get_processor_name(processor_name, &namelen);
  fprintf(stderr, "Process %d running on %s\n", myid, processor_name);

  string numbers = "0123456789";  // !!!!! np <= 10
  string myid_str(numbers, myid, 1);

  MPI_Status status;

  // define a new MPI data type for particles
  MPI_Datatype particletype;
  MPI_Type_contiguous(18, MPI_DOUBLE, &particletype);  // !!! 14->18 changed
  MPI_Type_commit(&particletype);

  //-------- end MPI init----------------

  // wait for gdb
  waitforgdb(myid);

  // read input file (e.g. patric.cfg):

  if(argv[1] == 0){
    printf("No input file name !\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  input_from_file(argv[1], myid);
  double eps_x = rms_emittance_x0;  // handy abbreviation
  double eps_y = rms_emittance_y0;  // same

  // Synchronous particle:

  SynParticle SP;
  SP.Z = Z;
  SP.A = A;
  SP.gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) ;
  SP.beta0 = sqrt((SP.gamma0*SP.gamma0-1.0)/(SP.gamma0*SP.gamma0)) ;
  SP.eta0 = 1.0/pow(gamma_t, 2)-1.0/pow(SP.gamma0, 2);

  //-------Init Lattice-------

  BeamLine lattice;
  double tunex, tuney;
  if(madx_input_file == 1)
    // read madx sectormap and twiss files
    lattice.init("../mad/", circum, tunex, tuney, SP.beta0, chroma);
  else{
    // init constant focusing (CF) sectormap and cell:
    SectorMap CF(CF_advance_h/NCF, CF_advance_v/NCF, CF_R, CF_length/NCF, SP.gamma0);
    BeamLine CF_cell;
    for(int j=0; j<NCF; j++)
      CF_cell.add_map(CF);
    lattice.init(CF_cell);
  }

  // Other variables:

  double dx = 2.0*piperadius/(NX-1.0);  // needed for Poisson solver and grids
  double dy = 2.0*piperadius/(NY-1.0);  // needed for Poisson solver and grids
  double dz = circum/NZ;
  double ds = 0.4;  // value needed here only for setting dxs, dys.
  double dxs = 4.0*(dx/ds)/(NX-1.0);  // only for plotting xs, not for tracking
  double dys = 4.0*(dx/ds)/(NX-1.0);  // only for plotting ys, not for tracking
  double charge = current*circum/(NPIC*SP.beta0*clight*qe);  // macro-particle charge Q/e
  double zm = 0.5*circum*bunchfactor;  // (initial) bunch length
  if(init_pic_z == 1 || init_pic_z == 3 || init_pic_z == 4 || init_pic_z == 6)
    zm = 1.5*0.5*circum*bunchfactor;  // for parabolic bunch
  double zm1 = -zm*1.0;  // left bunch boundary
  double zm2 = zm*1.0;  // right bunch boundary
  double rmsToFull;  // ratio of rms to full emittance for Bump; SP

  // open output file patric.dat:

  string data_dir = ausgabe;
  data_dir = data_dir + "/";
  string outfile = data_dir + "patric.dat";
  FILE *out = fopen(outfile.c_str(), "w");

  // init random number generator:

  long d = -11*(myid+1);  // was -1021  transverse distribution: each slice needs a different initialization !
  long dl = -103;  // was -103   longitudinal plane: same random set needed
  long dran = -101;  // for BTF noise excitation: same random sets needed


  // set some global lattice parameters

  double cell_length = lattice.get_L();
  int Nelements = lattice.get_size();
  if(myid == 0){
    cout << "Nelements:" << Nelements << endl;
    cout << "Cell length:" << cell_length << endl;
  }

  // define pointers to first/last element in beam line:

  const list<SectorMap>::iterator first_elem = lattice.get_first_element();
  const list<SectorMap>::iterator last_elem = --lattice.get_end_element();

  TwissP twiss0, twiss_TK;
  twiss0 = last_elem->get_twiss();
  twiss_TK = first_elem->get_twiss();
  double Ds0 = 0.0;  // Dispersion derivative

  if(madx_input_file == 0){
    // machine tunes from lattice
    lattice.phase_advance(tunex, tuney);
    tunex = circum/cell_length*tunex/(2.0*PI);
    tuney = circum/cell_length*tuney/(2.0*PI);
    if(myid == 0){
      cout << "advancex: " << tunex*180.0/PI << endl;
      cout << "tunex0: " << tunex << endl;
      cout << "tuney0: " << tuney << endl;
    }
  }

  // Chromatic correction kick:
  //  Chrom Chrom0(tunex, tuney, circum/(2.0*PI));

  // Octupole:
  Octupole Oct0(koct);

  // Amplitude detuning; works only for constant focusing; SP
  //  AmplitudeDetuning Amp0(tunex, tuney, dqx_detune/(1.0e-6*eps_x),
  //		 dqy_detune/(1.0e-6*eps_y), circum/(2.0*PI), CF);

  //--------end lattice----------

  // set matched RF voltage:

  int linrf = 0;
  if (cavity == 3) linrf = 1;
  double Ym = circum/(2.0*PI)*(1.0-cos(2.0*PI*zm/circum));
  if (linrf == 1) Ym = circum/(2.0*PI)*0.5*pow(2.0*PI*zm/circum, 2);
  double velm = abs(SP.eta0)*SP.beta0*clight*sqrt(5.0)*momentum_spread*
    2.0*PI/(circum);
  double fsyn = 1.0/(2.0*PI)*velm*sqrt(circum/(2.0*PI))/sqrt(2.0*Ym);
  double V0rf = pow(2.0*PI*fsyn, 2)*pow(circum, 2)/(2.0*PI)*mp*SP.A*
    SP.gamma0/(qe*SP.Z*abs(SP.eta0));
  if(myid == 0){
    cout << "Vrf [kV]: " << V0rf*1.0e-3 << "  fsyn [kHz]: " << fsyn*1.0e-3 << endl;
    // print IDL parameter file idl.dat:
    print_IDL(data_dir, numprocs, cell_length, Nelements);
  }

  // Init particle distribution:

  Pic Pics(&SP, charge, NPIC/numprocs, data_dir + "pics_" + myid_str + ".dat");
  Pics.z1 = zm1+myid*(zm2-zm1)/numprocs;  // left boundary in z for this slice
  Pics.z2 = Pics.z1+(zm2-zm1)/numprocs;  // right boundary
  double slice_length = Pics.z2-Pics.z1;  // slice length

  Pic NewPics(&SP, charge, NPIC/numprocs);
  NewPics.z1 = Pics.z1;
  NewPics.z2 = Pics.z2;

  // Init 1D longitudinal grids

  Grid1D rho_z_tmp(NZ, dz, -0.5*circum);
  Grid1D rho_z(NZ, dz, -0.5*circum, data_dir + "rho_z.dat");
  Grid1D dipole_current_x_tmp(NZ, dz, -0.5*circum);
  Grid1D dipole_current_x(NZ, dz, -0.5*circum, data_dir + "dipole_x.dat");
  Grid1D dipole_current_xs_tmp(NZ, dz, -0.5*circum);
  Grid1D dipole_current_xs(NZ, dz, -0.5*circum);
  Grid1D dipole_kick_x(NZ, dz, -0.5*circum, data_dir + "dipole_kick_x.dat");
  Grid1D dipole_current_y_tmp(NZ, dz, -0.5*circum);
  Grid1D dipole_current_y(NZ, dz, -0.5*circum, data_dir + "dipole_y.dat");

  // Init 2D transverse grids:

  Grid2D rho_xy(NX, NY, dx, dy, data_dir + "rho_xy.dat");
  Grid2D rho_xy_tmp(NX, NY, dx, dy);
  Grid2D xxs(NX, NX, dx, dxs, data_dir + "xxs.dat");
  Grid2D xxs_tmp(NX, NX, dx, dxs);
  Grid2D yys(NY, NY, dy, dys, data_dir + "yys.dat");
  Grid2D yys_tmp(NY, NY, dy, dys);
  Grid2D xsys(NX, NY, dxs, dys, data_dir + "xsys.dat");
  Grid2D xsys_tmp(NX, NY, dxs, dys);
  Grid2D zx(NZ, NX, dz, dx, data_dir + "zx.dat");
  Grid2D zx_tmp(NZ, NX, dz, dx);

  Grid2D Ex(NX, NY, dx, dy, data_dir + "Ex.dat");
  Grid2D Ey(NX, NY, dx, dy, data_dir + "Ey.dat");

  // Init 3D sliced grids (for 3D space charge calculation)

  if( fmod((float)NZ_bunch, (float)numprocs) != 0.0 ){
    cout << "NZ_bunch kein Vielfaches von numprocs" << endl;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  Grid3D rho_xyz(NZ_bunch/numprocs, Pics.z1, Pics.z2, rho_xy);
  Grid3D Ey3(NZ_bunch/numprocs, Pics.z1, Pics.z2, rho_xy);
  Grid3D Ex3(NZ_bunch/numprocs, Pics.z1, Pics.z2, rho_xy);

  // Init 2D Greens function for poisson solver

  Greenfb gf1(rho_xy, image_x, image_y);  // open boundary condition


  // set longitudinal distribution:

  // pointer to long. distribution function
  void (Pic::* long_dist)(double, double, double, long, long*, int);
  double d1 = 0., d2 = 0.;  // meaning depends on momentum distribution; SP
  int i1 = 0;  // same
  switch(init_pic_z){
  case 0:  //  coasting + Elliptic
    long_dist = &Pic::parabolic_dc;
    d1 = bunchfactor;
    d2 = circum;
    break;
  case 1:  //  bunch + Elliptic  (1.5 correction factor for bunching)
    long_dist = &Pic::parabolic;
    d1 = zm;
    break;
  case 2:  //  coasting + Gauss
    long_dist = &Pic::coast_gauss; //(bunchfactor, circum, momentum_spread, NPIC, &dl);
    d1 = bunchfactor;
    d2 = circum;
    break;
  case 3:  //  bunch + Gauss
    long_dist = &Pic::bunch_gauss;
    d1 = zm;
    d2 = circum;
    break;
  case 4:  //  const. bunch dist.
    long_dist = &Pic::bunch_const;
    d1 = zm;
    d2 = circum;
    i1 = linrf;
    break;
  case 5:  //  air bag dist.
    long_dist = &Pic::barrier_air_bag;
    d1 = zm;
    break;
  case 6:  //  bunch air bag dist.
    long_dist = &Pic::bunch_air_bag;
    d1 = zm;
    d2 = circum;
    break;
  default:
    printf("Invalid option for longitudinal particle distribution. Aborting.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }


  // set transverse distribution:

  // pointer to long. distribution function
  void (Pic::* trans_dist)(double, double, double, double, double, double, double,
			   double, double, double, long*);
  switch(init_pic_xy){
  case 0:  // Waterbag
    trans_dist = &Pic::waterbag_xy;
    rmsToFull = 6;
    break;
  case 1:  // KV
    trans_dist = &Pic::KV_xy;
    rmsToFull = 4;
    break;
  case 2:  // Semi-Gauss
    trans_dist = &Pic::SG;
    rmsToFull = 4;  // approximate
    break;
  case 3:  // Gauss
  trans_dist = &Pic::Gauss_xy;
  rmsToFull = 4;  // approximate
    break;
  default:
    printf("Invalid option for transverse particle distribution. Aborting.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  if(myid == 0)
    cout << "Expected single beamlett tune shifts: dQ_x="
	 << rp*SP.Z*current*circum / (rmsToFull*PI*clight*qe*SP.A*pow(SP.beta0*SP.gamma0, 3)*
				      (eps_x+sqrt(eps_x*eps_y*tunex/tuney)))*1e6
	 << ", dQ_y="
	 << rp*SP.Z*current*circum / (rmsToFull*PI*clight*qe*SP.A*pow(SP.beta0*SP.gamma0, 3)*
				      (eps_y+sqrt(eps_x*eps_y*tuney/tunex)))*1e6
	 << endl;


  // local orbit bump for beam injection; SP
  Bump lob(&lattice, eps_x, rmsToFull, momentum_spread, tunex, offcenter,
	   x_septum, inj_angle, max_inj, inj_phase, sept_return, myid);

  //----------------counters and other variables--------------------------

  int Nexchange = 1;  // exchange of particles between slices after every sector map.
  int Nprint = print_cell*Nelements;  // output of particles every cell*print_cell
  //int Nibs = 1;  // correct for IBS every Nibs steps
  double Ntot;  // total number of particles: for screen output
  int counter = 0;  // counts sector maps
  double s = 0.0;  // path length
  double Nslice;  // total number of slices
  double emitx;  // emittance: for screen output
  double dtheta = 0.0;  // btf dipole kick
  double pickup_h, pickup_v;  // horizontal/vertical pickup signals
  double rms_advancex = 0.0, rms_advancey = 0.0;  // rms phase advance: for output
  int inj_counter = 0;  // number of injected beamletts; SP
  long N_inj = 0;  // number of injected particles

  //---------parameters for exchange of particles between slices-------

  int destl;  //!< ID of left neighbour slice (-1: no neighbour).
  int destr;  //!< ID of right neighbour

  //---finite bunch: no exchange between ends---
  if(bc_end == 0){
    if(myid == 0){
      destl =-1;
      destr = myid+1;
    }else
      if(myid == numprocs-1){
	destl = myid-1;
	destr =-1;
      }else{
	destl = myid-1;
	destr = myid+1;
      }
  }
  //---periodic (in z) boundary condition---
  if(bc_end == 1){
    if(myid == 0){
      destl = numprocs-1;
      destr = myid+1;
    }else
      if(myid == numprocs-1){
	destl = myid-1;
	destr = 0;
      }
      else{
	destl = myid-1;
	destr = myid+1;
      }
  }

  //--------------------- end-parameters for particle exchange ---------------

  long *septLoss = new long;
  long *sl_slice = new long;
  double *momenta = new double[17];
  double *momenta_tot = new double[17];

  //--------------------------------------------------------------------------
  //----------------------- start loop (do...while) --------------------------
  //--------------------------------------------------------------------------

  do{  // injection; SP
    if(!(counter%Nelements)){  // at beginning each turn...
      if(inj_counter < max_inj){  // ... inject beamlett (if apropriate) ...
	(NewPics.*long_dist)(d1, d2, momentum_spread, NPIC, &dl, i1);
	(NewPics.*trans_dist)(1.e-6*eps_x, 1.0e-6*eps_y, twiss_TK.alpx, twiss_TK.alpy,
			      pow(mismatch_x, 2)*twiss_TK.betx, pow(mismatch_y, 2)*twiss_TK.bety,
			      twiss_TK.Dx, Ds0, offcenter, inj_angle, &d);
	*sl_slice = NewPics.localLoss_x(x_septum, 100.);  // loss on septum
	MPI_Reduce(sl_slice, septLoss, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	++inj_counter;
	Pics.add_particles(*NewPics.get_particles());
	NewPics.clear_particles();
	N_inj += NPIC;
	if(myid == 0)
	  cout<<"The incoming beamlett number "<<inj_counter<< " lost "<<*septLoss<<
	    " macro particles on the septum.\n";
      }
      lob.decrement();  // ... reduce height of orbit bump
    }

    //------------ Start Output----------------------------------------
    
    // store rms momenta every time step in patric.dat:
    
    if(counter%1 == 0){
      Nslice = Pics.get_size();  // number of particles in this slice
      momenta[0] = Nslice*Pics.rms_emittance_x();
      momenta[1] = Nslice*Pics.rms_emittance_y();
      momenta[2] = Nslice*Pics.x_max();
      momenta[3] = Nslice*Pics.y_max();
      momenta[4] = Nslice*Pics.x_rms();
      momenta[5] = Nslice*Pics.y_rms();
      momenta[6] = Nslice*Pics.rms_momentum_spread();
      momenta[7] = Nslice*Pics.xzn(2.0, zm);
      momenta[8] = Nslice*Pics.xzn(1.0, zm);
      momenta[9] = Nslice;
      momenta[10] = Nslice*rms_advancex;  // rms phase advance in x
      momenta[11] = Nslice*rms_advancey;
      momenta[12] = Nslice*Pics.offset_x();
      momenta[13] = Nslice*Pics.offset_y();
      momenta[14] = Nslice*dtheta;  // btf noise signal
      momenta[15] = Nslice*pickup_h;
      momenta[16] = Nslice*pickup_v;
      // mpi_reduce for summation of all 17 moments over all slices
      MPI_Reduce(momenta, momenta_tot, 17, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      Ntot = momenta_tot[9];  // total number of particles over all slices
      emitx = momenta_tot[0]/Ntot;  // total rms emittance

      // stop when loss tolerance level is exceeded
      if(myid == 0 && Ntot/N_inj <= lossTol){  // test on numer of injected particles; SP
	cout<<"Loss tolerance exceeded within "<<counter/Nelements+1<<" turns ("<<
	  Ntot<<" of "<<N_inj<<" macro particles left). Exiting.\n";
	cout.flush();
	MPI_Abort(MPI_COMM_WORLD, 0);
      }
      //      cout<<counter<<' '<<lattice.get_element()->get_name()<<' '<<lattice.get_element()->get_K(1)<<endl;  //tmp
      // write momenta
      if(myid == 0){
	fprintf(out, "%g", s);
	for(int i=0; i<17; i++)
	  if(i != 9)
	    fprintf(out, "%15g", momenta_tot[i]/Ntot);
	  else
	    fprintf(out, "%15g", momenta_tot[i]);
	fprintf(out, "\n");
	fflush(out);
      }
    }


    //------output every Nprint*sectormap---------

    if(counter%Nprint == 0){
      if(myid == 0){
	// to screen
	printf("saving at s=%g (m) eps_t=%g dp/p=%g zm2=%g Ntotal=%g\n", s, 1.0e6*emitx,
	       Pics.rms_momentum_spread(), zm2, Ntot);
	cout.flush();
	
	// electric fields
	Ex.print();
	Ey.print();
      }

      // paricle coordinates to pic.dat:
      Pics.print(pic_subset);

      // collect densities for output only:

      Pics.gatherZ(charge*qe/dz, rho_z_tmp);
      Pics.gatherX(SP.beta0*clight*charge*qe/dz, dipole_current_x_tmp);
      Pics.gatherY(SP.beta0*clight*charge*qe/dz, dipole_current_y_tmp);
      Pics.gatherXY(charge*qe/circum, rho_xy_tmp);
      Pics.gatherXXs(charge*qe/circum, xxs_tmp);
      Pics.gatherYYs(charge*qe/circum, yys_tmp);
      Pics.gatherXsYs(charge*qe/circum, xsys_tmp);
      Pics.gatherZX(charge*qe/circum, zx_tmp);

      // summation over all slices:

      MPI_Allreduce(rho_z_tmp.get_grid(), rho_z.get_grid(),
		    NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(dipole_current_x_tmp.get_grid(), dipole_current_x.get_grid(),
		    NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(dipole_current_y_tmp.get_grid(), dipole_current_y.get_grid(),
		    NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(rho_xy_tmp.get_grid(), rho_xy.get_grid(),
		    NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(xxs_tmp.get_grid(), xxs.get_grid(),
		    NX*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(yys_tmp.get_grid(), yys.get_grid(),
		    NY*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(xsys_tmp.get_grid(), xsys.get_grid(),
		    NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(zx_tmp.get_grid(), zx.get_grid(),
		    NZ*NX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
      // output to density files:

      if(myid == 0){
	dipole_current_x.print();
	dipole_kick_x.print();
	dipole_current_y.print();
	rho_z.print();
	rho_xy.print();
	xxs.print();
	yys.print();
	xsys.print();
	zx.print();
      }
    }

    //-----------------end output--------------------------------------------


    // at beginning of a cell: calculate advance per (last) cell,
    // store old coordinates

    if(lattice.get_element() == first_elem){
      rms_advancex = Pics.rms_phaseadvance_h();  // Pics.rms_wavelength_h();
      rms_advancey = Pics.rms_phaseadvance_v();  // Pics.rms_wavelength_v();
      if(footprint == 0)
	Pics.store_old_coordinates();
    }

    // Transport particles through sectormap, update slice position s:

    ds = lattice.get_element()->get_L();
    s += ds;
    if(lattice.get_element()->get_name() == "\"SEPTUM\"")  // losses at septum; SP
      Pics.localLoss_x(-piperadius, coll_halfgap);

    double xi_h = lattice.get_element()->get_twiss().xix;
    double xi_v = lattice.get_element()->get_twiss().xiy;
    Pics.transport(lattice.get_element()->get_map(), xi_h, xi_v, piperadius);

    //-----exchange particles between slices------------------------

    if(counter != 0 && counter%Nexchange == 0 && numprocs > 1){
      int Npl;  //!< Number of particles to be exchanged with left neighbour
      int Npr;  //!< particles exchanged with right neighbour

      //! vector of particles to be exchanged
      vector<Particle> pl, pr;

      // send particle to neighbor slices:
	
      if(destl >= 0){
	pl = Pics.get_particles_left(circum);
	Npl = pl.size();
	MPI_Send(&Npl, 1, MPI_INT, destl, 1, MPI_COMM_WORLD);
	MPI_Send(&pl[0], Npl, particletype, destl, 1, MPI_COMM_WORLD);
      }
      if(destr >= 0){
	pr = Pics.get_particles_right(circum);
	Npr = pr.size();
	MPI_Send(&Npr, 1, MPI_INT, destr, 0, MPI_COMM_WORLD);
	MPI_Send(&pr[0], Npr, particletype, destr, 0, MPI_COMM_WORLD);
      }

      // receive from neighbour slices:
	
      Npl = 0; Npr = 0;
      vector<Particle> pl_in, pr_in;
      if( destl >= 0 ){
	MPI_Recv(&Npl, 1, MPI_INT, destl, 0, MPI_COMM_WORLD, &status);
	pl_in = vector<Particle>(Npl);
	MPI_Recv(&pl_in[0], Npl, particletype, destl, 0, MPI_COMM_WORLD, &status);
      }
      if(destr >= 0){
	MPI_Recv(&Npr, 1, MPI_INT, destr, 1, MPI_COMM_WORLD, &status);
	pr_in = vector<Particle>(Npr);
	MPI_Recv(&pr_in[0], Npr, particletype, destr, 1, MPI_COMM_WORLD, &status);
      }
      Pics.add_particles(pl_in);
      Pics.add_particles(pr_in);
    }

    //-----end exchange of particles-------------

    // periodic bc without exchange
    if(numprocs == 1)
      Pics.periodic_bc(circum);	

    // update wave lengths

    //if( footprint == 1){
    //Pics.update_wavelength_h(ds, 0.0);
    //Pics.update_wavelength_v(ds);}

    // nonlinear thin lens kick:
    if(octupole_kick == 1)
      Pics.kick(Oct0, ds);

    //if(ampdetun_kick == 1)  // works only for constant focusing
    //Pics.kick(Amp0, ds);

    // correct for chromaticity
    //    if(chroma == 1)  // chromaticity for CF disabled; SP
    //      Pics.kick(Chrom0, ds);

    // cavity kick every cell:

    if(cavity == 1 && counter%Nelements == 0.0)
      Pics.cavity_kick(V0rf*cell_length/circum, 1, circum/(2.0*PI));
    if(cavity == 2 && counter%Nelements == 0.0)
      Pics.barrier_kick(zm1, zm2);
    if(cavity == 3 && counter%Nelements == 0.0)
      Pics.cavity_kick_linear(V0rf*cell_length/circum, 1, circum/(2.0*PI));

    // Pickup signals

    Pics.gatherX(SP.beta0*clight*charge*qe/dz, dipole_current_x_tmp);
    Pics.gatherY(SP.beta0*clight*charge*qe/dz, dipole_current_y_tmp);
    MPI_Allreduce(dipole_current_x_tmp.get_grid(), dipole_current_x.get_grid(), NZ,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
    MPI_Allreduce(dipole_current_y_tmp.get_grid(), dipole_current_y.get_grid(), NZ,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    pickup_h = Pics.pickup_signal(dipole_current_x, circum,
				  s/(SP.beta0*clight))/current;
    pickup_v = Pics.pickup_signal(dipole_current_y, circum,
				  s/(SP.beta0*clight))/current;


    //---------------impedance kicks-----------------------

    komplex dqc_t(dqcr, dqci);  // for sliced == 0

    if(imp_kick == 1){
      if(sliced == 0)
	Pics.kick(ds/circum*InducedKick(Pics.offset_x(), ds, dqc_t, SP.beta0,
					tunex, circum), 0.0);
      else{
	dipole_kick_x.reset(); 	
	if(Rs > 0.0 || leit > 0.0){
	  Pics.gatherXs(SP.beta0*clight*charge*qe/dz, dipole_current_xs_tmp);
	  MPI_Allreduce(dipole_current_xs_tmp.get_grid(), dipole_current_xs.get_grid(),
			NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
	  InducedWakeKick(dipole_kick_x, dipole_current_x, dipole_current_xs, tunex,
			  2.0*PI*SP.beta0*clight/circum, nres, Rs, Qs, piperadius,
			  leit, SP.beta0, SP.gamma0*mp*SP.A*pow(clight, 2), SP.Z*qe);
	}
	if(Zimage != 0.0)
	  InducedKick(dipole_kick_x, dipole_current_x, Zimage, SP.beta0,
		      SP.gamma0*mp*SP.A*pow(clight, 2), SP.Z*qe);
	Pics.impedance_kick(dipole_kick_x, circum, ds);
      }
    }

    //---------------end impedance kicks-----------------------

    //------------self-consistent space charge kicks after every sectormap----

    if(space_charge == 1){
      // PIC -> charge density for Poisson solver:
      if (sliced == 0){
	Pics.gatherXY(charge*qe/circum, rho_xy_tmp);		
	MPI_Allreduce(rho_xy_tmp.get_grid(), rho_xy.get_grid(), NX*NY, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);	
      }else{	
	Pics.gatherXYZ(charge*qe/rho_xyz.get_dz(), rho_xyz);

	// send and receive density ghost grids to neighbor slices:
	// what is exchanged here ???
	if(destl >= 0)
	  MPI_Send(rho_xyz.get_ghostl(), NX*NY, MPI_DOUBLE, destl, 2, MPI_COMM_WORLD);
	if(destr >= 0){
	  MPI_Recv(rho_xy_tmp.get_grid(), NX*NY, MPI_DOUBLE, destr, 2, MPI_COMM_WORLD,
		   &status);
	  rho_xyz[NZ_bunch/numprocs-1] += rho_xy_tmp;
	}
	if(destr >= 0)
	  MPI_Send(rho_xyz.get_ghostr(), NX*NY, MPI_DOUBLE, destr, 3, MPI_COMM_WORLD);
	if(destl >= 0){
	  MPI_Recv(rho_xy_tmp.get_grid(), NX*NY, MPI_DOUBLE, destl, 3, MPI_COMM_WORLD,
		   &status);
	  rho_xyz[0]+= rho_xy_tmp;
	}
      }
	
      // Poisson solver
      if(sliced == 0)
	poisson_xy(Ex, Ey, rho_xy, gf1);
      else{
	poisson_xyz(Ex3, Ey3, rho_xyz, gf1);
	
	// send and receive efield ghost grids to neighbor slices:
	if(destl >= 0){
	  MPI_Send(Ex3.get_ghostl(), NX*NY, MPI_DOUBLE, destl, 2,
		   MPI_COMM_WORLD);
	  MPI_Send(Ey3.get_ghostl(), NX*NY, MPI_DOUBLE, destl, 4,
		   MPI_COMM_WORLD);
	}
	if(destr >= 0){
	  MPI_Recv(Ex3[NZ_bunch/numprocs-1].get_grid(), NX*NY, MPI_DOUBLE,
		   destr, 2, MPI_COMM_WORLD, &status);
	  MPI_Recv(Ey3[NZ_bunch/numprocs-1].get_grid(), NX*NY, MPI_DOUBLE,
		   destr, 4, MPI_COMM_WORLD, &status);
	}
	if(destr >= 0){
	  MPI_Send(Ex3.get_ghostr(), NX*NY, MPI_DOUBLE, destr, 3, MPI_COMM_WORLD);
	  MPI_Send(Ey3.get_ghostr(), NX*NY, MPI_DOUBLE, destr, 5, MPI_COMM_WORLD);
	}
	if(destl >= 0){
	  MPI_Recv(Ex3[0].get_grid(), NX*NY, MPI_DOUBLE, destl, 3, MPI_COMM_WORLD, &status);
	  MPI_Recv(Ey3[0].get_grid(), NX*NY, MPI_DOUBLE, destl, 5, MPI_COMM_WORLD, &status);
	}
      }
    }
    
    // Shift xs and ys:
    
    if(space_charge == 1 && ds > 0.0){
      if(sliced == 0)
	Pics.kick(Ex, Ey, ds);
      else
	Pics.kick(Ex3, Ey3, ds);
    }
    
    //---------------end self-consistent space charge kicks---------------

    // linear sc kicks:
    
    if(space_charge == 2 && ds > 0.0)
      Pics.linear_SC_kick(dQxm, dQym, tunex, tuney, rho_z, current/(SP.beta0*clight),
			  dipole_current_x, dipole_current_y, circum, ds);
	
    // nonlinear sc kicks:

    if(space_charge == 3 && ds > 0.0)
      Pics.nonlinear_SC_kick(sqrt(1.0e-6*twiss0.betx*eps_x), sqrt(1.0e-6*twiss0.bety*eps_y),
			     dQxm, dQym, tunex, tuney, rho_z, current/(SP.beta0*clight),
			     circum, ds);

    // dipole noise modulation kick:
    double dnoiseamp = 1.0e-6;
    double nus = fsyn/(SP.beta0*clight/circum);
    if(btf == 1)
      dtheta = Pics.dipole_mod_kick(s/(SP.beta0*clight), ds, circum, dnoiseamp,
				    (tunex+nus)*SP.beta0*clight/circum, btf_harmonic);	
						
    // correct for ibs:

    /*if(counter != 0 && counter%Nibs == 0){
      double rate_ibs = 1.0e4;
      double Dz = rate_ibs*pow(Pics.rms_momentum_spread(), 2);
      double Dxy = rate_ibs*0.5*(Pics.rms_emittance_x()+Pics.rms_emittance_y());
      double betx = lattice.get_element()->get_betx();
      double bety = lattice.get_element()->get_bety();
      Pics.langevin(rate_ibs, rate_ibs*0.0, Dxy, Dz*0.0, Nibs*ds, betx, bety,
        &d);
      }*/

    // For bunch compression: Update slice boundaries z1 and z2 from
    // new bunch boundaries zm1, zm2:

    /*if(counter != 0 && counter%Nexchange == 0){
      if(myid == 0)
      zm1 = Pics.z_min();
      MPI_Bcast(&zm1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(myid == numprocs-1)
      zm2 = Pics.z_max();
      MPI_Bcast(&zm2, 1, MPI_DOUBLE, numprocs-1, MPI_COMM_WORLD);

      Pics.z1 = zm1+myid*(zm2-zm1)/numprocs;
      Pics.z2 = Pics.z1+(zm2-zm1)/numprocs;
      slice_length = Pics.z2-Pics.z1;

      rho_xyz.get_zleft() = zm1;
      rho_xyz.get_zright() = zm2;
      Ex3.get_zleft() = zm1;
      Ex3.get_zright() = zm2;
      Ey3.get_zleft() = zm1;
      Ey3.get_zright() = zm2;
      }*/


    // advance in beam line, go to next element:

    lattice.next_element();
    ++counter;

  }while(counter != cells*Nelements);

  //------------------end of loop-------------------------------

  // close files, free heap:

  delete septLoss, sl_slice;
  delete[] momenta, momenta_tot;  // [] needed here!; SP
  fclose(out);

  // MPI end:

  MPI_Finalize();

  time2 = time(0);
  double sec = difftime(time2, time1);
  double h = floor(sec/3600);
  double min = floor(sec/60-60.*h);
  sec -= 3600.*h+60.*min;
  if(myid == 0)
    cout << "Total losses: " << (1-Ntot/(max_inj*NPIC))*100. << " \%\n" <<
      "Stored particles: " << current*circum*Ntot/(qe*Z*SP.beta0*clight*NPIC) << endl <<
      "Computation time: " << h << ":" << min << ":" << sec << endl;
}
