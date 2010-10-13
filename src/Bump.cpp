// This class defines a local orbit bump for the injection of the beam.
// 4 kickers and a marker called INJ for the injection must be contained
// in the madx twiss file.
// created by Stefan Paret (SP), May 2010

using namespace std;

#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <mpi.h>
#include <PhyConstants.h>

typedef vector<double> vektor;

#include "Bump.h"
#include "SectorMap.h"


Bump::Bump(BeamLine* bl, double emit_x, double rmsToFull, double dp0, double Q_x,
	   double &offset, double x_septum, double inj_angle, unsigned &max_inj,
	   double inj_phase, unsigned sept_return, int id){
  double d_septum = 0.0001;  // septum thickness
  list<SectorMap>::iterator el = bl->get_first_element();
  const list<SectorMap>::iterator ending = bl->get_end_element();

  // adjust injection settings to beam parameters
  double a = sqrt(el->get_betx()*emit_x*rmsToFull)*0.001+el->get_twiss().Dx*dp0;  // half width of injected beam [m] with WB distribution
  offset += x_septum + d_septum + a;  // offset relative to beam center -> w.r.t. nominal orbit

  // find injection point, injection kickers and their lattice functions
  double beta_I, psi_I, alpha_I, psi[4], beta[4], defl[4];
  while(el != ending){
    if(el->get_name() == "\"INJ\""){
      beta_I = el->get_betx();
      psi_I = el->get_mux();
      alpha_I = el->get_alpx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S01MB3\""){
      kickers[0] = &*el;
      psi[0] = el->get_mux()*2*PI;
      beta[0] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S03MB4\""){
      kickers[1] = &*el;
      psi[1] = el->get_mux()*2*PI;
      beta[1] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S11MB1\""){
      kickers[2] = &*el;
      psi[2] = (el->get_mux()-Q_x)*2*PI;
      beta[2] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S12MB2\""){
      kickers[3] = &*el;
      psi[3] = (el->get_mux()-Q_x)*2*PI;
      beta[3] = el->get_betx();
      break;
    }
    ++el;
  }
  if(el == ending){
    cout<<"Error: Not all kickers for local orbit bump could be found or mark for injection point is missing.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  // The bump height is defined by user given offcenter parameter. The injection angle is equal to the septum tilt angle (as done in SIS18).
  const double amp0 = offset;
  const double ampp0 = inj_angle;
  const double delAmp = (amp0-2*a)/double(max_inj);

  // calculate deflection angles for local orbit bump adjusted to incoming beam
  double b_21 = sqrt(beta[3]*beta[2])*sin(psi[3]-psi[2]);
  double b_I1 = sqrt(beta_I*beta[2])*sin(psi_I-psi[2]);
  double b_I2 = sqrt(beta_I*beta[3])*sin(psi_I-psi[3]);
  double b_31 = sqrt(beta[0]*beta[2])*sin(psi[0]-psi[2]);
  double b_32 = sqrt(beta[0]*beta[3])*sin(psi[0]-psi[3]);
  double b_41 = sqrt(beta[1]*beta[2])*sin(psi[1]-psi[2]);
  double b_42 = sqrt(beta[1]*beta[3])*sin(psi[1]-psi[3]);
  double b_43 = sqrt(beta[1]*beta[0])*sin(psi[1]-psi[0]);
  double d_I1 = sqrt(beta[2]/beta_I)*cos(psi_I-psi[2]) - alpha_I/beta_I*b_I1;
  double d_I2 = sqrt(beta[3]/beta_I)*cos(psi_I-psi[3]) - alpha_I/beta_I*b_I2;

  defl[2] = (d_I2*amp0 - b_I2*ampp0)/b_21;  // S11MB1
  defl[3] = (-d_I1*amp0 + b_I1*ampp0)/b_21;  // S12MB2
  defl[0] = (-b_41*defl[2] - b_42*defl[3])/b_43;  // S01MB3
  defl[1] = (b_31*defl[2] + b_32*defl[3])/b_43;  // S03MB4

  if(id == 0){
    cout<<"Injection parameters:"<<
      "a/m="<< a <<", offset/m=" << offset <<", amp0/m=" << amp0 << ", delAmp/m=" << delAmp <<
      ", ampp0= " << ampp0<< ", Q_x=" << Q_x << ", max_inj=" << max_inj << endl;
    //cout << "beta_I/m=" << beta_I << ", alpha_I=" << alpha_I << endl;
    cout << "Deflection angles/(mm mrad): ";
    for(short i=0; i<4; ++i)
      cout << defl[i]*1000. << ", ";
    cout << endl;
    for(short i=0; i<4; ++i){  // check against maximal deflection
      if(defl[i] > 0.0084)
	cout << "WARNING: Required deflection angle in kicker " << i
	     << " exceeds limit of 8.4 mrad.\n";
    }
    if(max_inj<1){
      cout<<"Error: Beam too large for injection. Execution aborted.\n";
      cout.flush();
      MPI_Abort(MPI_COMM_WORLD, 0);
    }
  }

  double ab = delAmp/amp0;  // relative decrement
  for(short i=0; i<4; ++i){
    decr[i] = ab*defl[i];  // set reduction per revolution
    kickers[i]->get_K(1) = defl[i];  // set deflection angle
  }

  // compensate first decrementation at beginning of loop in Main
  kickers[0]->get_K(1) += decr[0];
  kickers[1]->get_K(1) += decr[1];
}


void Bump::decrement(){
  // works only as long as beam line in Main is not moved
  if(kickers[0]->get_K(1) != 0.)
    if(kickers[0]->get_K(1) > decr[0])
      for(short i=0; i<4; ++i)
	kickers[i]->get_K(1) -= decr[i];
    else
      for(short i=0; i<4; ++i)
	kickers[i]->get_K(1) = 0.;
}
