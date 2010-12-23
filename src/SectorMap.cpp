#include <string>
#include <list>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

typedef vector<double> vektor;

#include "PhyConstants.h"
#include "SectorMap.h"

/*!
  Generates a linear transport map for a simple continous focusing element
  of length = length, bending radius R and phase advances [grad] sigma_x and sigma_y
*/

SectorMap::SectorMap(double sigx, double sigy, double R, double length,
		     double gamma0){
  ElementName = "CF";
  L = length;
  double beta0 = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0));
  double eta0 = pow(R*sigx/length, 2)-1.0/pow(gamma0, 2);
  for(int j = 0; j<36; j++) T[j] = 0.0;
  for(int j = 0; j<6; j++) K[j] = 0.0;
  get_T(0, 0) = cos(sigx);
  get_T(0, 1) = length/sigx*sin(sigx);
  if (R>0.0)
    get_T(0, 5) = length*length/(sigx*sigx*R)*(1.0-cos(sigx))/beta0;
  else
    get_T(0, 5) = 0.0;
  get_T(1, 0) =-sigx/length*sin(sigx);
  get_T(1, 1) = cos(sigx);
  if (R>0.0)
    get_T(1, 5) = length/(sigx*R)*sin(sigx)/beta0;
  else
    get_T(1, 5) = 0.0;
  get_T(2, 2) = cos(sigy);
  get_T(2, 3) = length/sigy*sin(sigy);
  get_T(3, 2) =-sigy/length*sin(sigy);
  get_T(3, 3) = cos(sigy);
  get_T(4, 4) = 1.0;
  get_T(4, 5) =-eta0*length/(beta0*beta0);
  get_T(5, 5) = 1.0;
  twiss.betx = length/sigx;
  twiss.bety = length/sigy;
  twiss.alpx = 0.0;
  twiss.alpy = 0.0;
  if(R > 0.0)
    twiss.Dx = pow(length/sigx, 2)/R;
  else
    twiss.Dx = 0.0;
}


SectorMap& SectorMap::operator = (const SectorMap& M){
  if( this != &M ){
    ElementName = M.ElementName;	
    L = M.L;
    for(int j = 0; j<36; j++) T[j] = M.T[j];
    for(int j = 0; j<6; j++) K[j] = M.K[j];
    twiss = M.twiss;
  }
  return *this;
}
	

SectorMap SectorMap::operator*(const SectorMap& M){
  SectorMap tmap;
  for(int i=0; i<6; ++i)
    for(int j=0; j<6; ++j)
      for(int l=0; l<6; ++l)
	tmap.T[j*6+i] += T[j*6+l]*M.T[l*6+i];
  tmap.L = L+M.L;
  return tmap;
}


/*!
Detailed description
*/

void SectorMap::transport(vektor& R1, vektor& R0, double K_xix, double K_xiy){
  for(int j=0; j<6; ++j){
    R1[j] = K[j];  // was 0.0, K[j] inserted for bumper kicks; SP
    for(int l=0; l<6; ++l)
      R1[j] += T[j*6+l]*R0[l];
  }
  R1[1] += K_xix*R0[0]*R0[5];  // chromaticity kick; SP
  R1[3] += K_xiy*R0[2]*R0[5];  // chromaticity kick; SP
}


/*!
  Calculates the phase advances per element.
*/

void SectorMap::phase_advance(double& sigx, double& sigy){
  double ah = get_T(0, 0);
  double dh = get_T(1, 1);
  double av = get_T(2, 2);	
  double dv = get_T(3, 3);
  sigx = acos(0.5*(ah+dh));
  sigy = acos(0.5*(av+dv));
}


/*!
  Constructs a beam line from MADX twiss and sectormap file in the
  dir directory. Sets the iterator to the first element.
*/

void BeamLine::init(string dir, double &length, double &Q_hor, double &Q_ver, double beta, int chroma){
  read_madx_twiss(dir + "twiss_inj.txt", length, Q_hor, Q_ver, beta, chroma);
  read_madx_sectormap(dir + "sectormap_inj.txt");
  element = line.begin(); 		
}


/*!
  Reads the MADX twiss file and generates the beam line.
  This function needs to be called first.
*/

void BeamLine::read_madx_twiss(string fname, double &circum, double &Q_hor, double &Q_ver,
			       double beta, int chroma){
  int index;
  string str;
  SectorMap SMap;
  double s, l;
  TwissP tw;
  list<string> pars;
  double ax0, ay0, bx0, by0, l0, mux0, muy0, tmp1;

  ifstream twissfile(fname.c_str());
  if(twissfile.fail()){
    cout<<"Twiss file could not be opened.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  do{  // read header
    getline(twissfile, str);
    if(str.find("@ LENGTH") != string::npos){
      index = str.find_first_of("0123456789");
      circum = atof(str.substr(index).c_str());  // keep beam length
    }
    if(str.find("@ Q1") != string::npos){
      index = str.find_first_of("0123456789", 5);
      Q_hor = atof(str.substr(index).c_str());  // keep horizontal tune
    }
    if(str.find("@ Q2") != string::npos){
      index = str.find_first_of("0123456789", 5);
      Q_ver = atof(str.substr(index).c_str());  // keep horizontal tune
    }
  }while( str.find("*") == -1 );
  int first = 50, last;  // identify parameters in file
  while( (first = str.find_first_not_of(' ', first)) != string::npos ){
    last = str.find_first_of(' ', first);
    pars.push_back(str.substr(first, last-first));
    first = last+1;
  }

  // adjust pointers to point on column with data corresponding to name
  int numPars = pars.size();
  double **target = new double*[numPars];
  list<string>::iterator parsIter = pars.begin();
  for(int i=0; i<numPars; ++i){
    if(*parsIter=="S")
      target[i] = &s;
    else if(*parsIter=="L")
      target[i] = &l;
    else if(*parsIter=="ALFX")
      target[i] = &tw.alpx;
    else if(*parsIter=="ALFY")
      target[i] = &tw.alpy;
    else if(*parsIter=="BETX")
      target[i] = &tw.betx;
    else if(*parsIter=="BETY")
      target[i] = &tw.bety;
    else if(*parsIter=="DX")
      target[i] = &tw.Dx;
    else if(*parsIter=="MUX")
      target[i] = &tw.mux;
    else if(*parsIter=="MUY")
      target[i] = &tw.muy;
    else if(*parsIter=="DMUX")
      target[i] = &tw.xix;
    else if(*parsIter=="DMUY")
      target[i] = &tw.xiy;
    else
      target[i] = 0;
    ++parsIter;
  }

  getline(twissfile, str);  // skip line with data format information
  getline(twissfile, str);  // skip first line with numbers: repeated at end of file
  do{
    twissfile>>SMap.get_name()>>str;  // read strings in first 2 columns
    // first two columns have to be name and keyword, otherwise reading fails
    for(int i=0; i<numPars; ++i)  // read parameters
      twissfile>>*target[i];
    SMap.get_L() = l;
    // postprocess input
    tw.Dx *= beta;  // MAD's DX -> D_x
    tw.mux *= 2*PI;  // tune to phase
    tw.muy *= 2*PI;  // tune to phase
    if(chroma && str == "\"QUADRUPOLE\""){  // if chromaticity is on and beam line element is quadrupole set constants for chromaticity kicks
      // for explanation of chromaticity kicks see notes p. 148
      tmp1 = sqrt(tw.betx/bx0)*(cos(tw.mux-mux0)+ax0*sin(tw.mux-mux0));
      if(tmp1 < 1.)
	tw.xix = pow(acos(tmp1), 2)/l/beta;  // k_quad*l_quad; x focusing
      else
	tw.xix = -pow(acosh(tmp1), 2)/l/beta;  // k_quad*l_quad; x defocusing
      tw.xiy = -tw.xix;  // opposite behaviour in y
    }else{  // no chromaticity or not in quadrupole -> no kick
      tw.xix = 0.;
      tw.xiy = 0.;
    }
    ax0 = tw.alpx;  // keep parameters for (possible) chromaticity kick
    bx0 = tw.betx;
    mux0 = tw.mux;
    SMap.get_twiss() = tw;
    line.push_back(SMap);
  }while( SMap.get_name().find("$END") == -1 );
  if(twissfile.fail()){
    cout<<"Error reading twiss file.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  twissfile.close();
  delete[] target;
}


/*!
  Initializes all sectormaps in the beam line.
  This function must be called after read_madx_twiss.
*/

void BeamLine::read_madx_sectormap(string fname){
  string str;	
  int j, l, u;
  list<SectorMap>::iterator pos = line.begin();
  double ddummy;
  string sdummy;

  ifstream mapfile(fname.c_str());
  if(mapfile.fail()){
    cout<<"Sectormap file could not be opened.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  do
    getline(mapfile, str);
  while( str.find("$START") == -1 );

  for(u = 0; u < line.size(); ++u){
    mapfile >> sdummy >> ddummy;
    for(j = 0; j < 6; ++j)
      mapfile >> pos->get_K(j);
    for(l = 0; l < 6; ++l)
      for(j = 0; j < 6; ++j)
	mapfile >> pos->get_T(j, l);
    for(l = 0; l < 36; ++l)
      for(j = 0; j < 6 ; ++j)
	mapfile >> ddummy;
    ++pos;
  }

  if(mapfile.fail()){
    cout<<"An error occured when reading sectormap file.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  mapfile.close();
}


double BeamLine::get_L(){
  double tem = 0.0;
  list<SectorMap>::iterator tpos = line.begin();
  for(int j = 0; j< line.size(); ++j){
    tem += tpos->get_L();
    ++tpos;
  }
  return tem;
}
	

void BeamLine::next_element(){
  if(element != line.end())
    ++element;
  if(element == line.end())
    element = line.begin();
}


void  BeamLine::set_element(int j){
  if(j < 0 || j >= line.size()){
    cout << "BeamLine::set_element: Error !" << endl;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  element = line.begin();
  for(int i = 0; i<j; i++)
    next_element();
}


/*!
  Returns a long linear sectormap by
  multiplication of the sectormaps
  in the beam line from element start
  to element end.
*/

/*
  SectorMap BeamLine::submap(int start, int end){
  set_element(start);
  SectorMap tmap(*element);
  SectorMap tmap2;
  next_element();
  for(int u = start+1; u < end; u++){
  tmap2 = element->operator*(tmap);
  tmap = tmap2;
  if ( u == end-1 ) tmap.get_twiss() = element->get_twiss();
  next_element();
  }
  return tmap;
  }
*/

/*!
  Phase advances calculated from elements.
*/


void BeamLine::phase_advance(double& sigx, double& sigy){
  list<SectorMap>::iterator pos = line.begin();
  SectorMap tmap(*pos->get_map());
  SectorMap tmap2;
  ++pos;
  list<SectorMap>::iterator pos0 = pos;
  for(pos = pos0; pos != line.end(); pos++){
    tmap2 = pos->operator*(tmap);
    tmap = tmap2;
  }	
  sigx = acos(0.5*(tmap.get_T(0, 0)+tmap.get_T(1, 1)));
  sigy = acos(0.5*(tmap.get_T(2, 2)+tmap.get_T(3, 3)));
}


void BeamLine::print_twiss(){
  list<SectorMap>::iterator pos = line.begin();
  list<SectorMap>::iterator pos0 = pos;
  SectorMap tmap;
  TwissP ttw;
  for(pos = pos0; pos != line.end(); pos++){
    ttw = pos->get_twiss();
    cout << "name:" << pos->get_name() << endl;
    cout << "betax:" << ttw.betx << endl;
    cout << endl << flush;
  }
}


void BeamLine::print_T(){
  list<SectorMap>::iterator pos = line.begin();
  list<SectorMap>::iterator pos0 = pos;
  double tmap;
  for(pos = pos0; pos != line.end(); pos++){
    cout << "name:" << pos->get_name() << endl;
    cout << pos->get_T(0, 0) << " " << pos->get_T(0, 1) << " " << pos->get_T(0, 2) << endl;
    cout << pos->get_T(1, 0) << " " << pos->get_T(1, 1) << " " << pos->get_T(1, 2) << endl;
    cout << pos->get_T(2, 0) << " " << pos->get_T(2, 1) << " " << pos->get_T(2, 2) << endl;
    cout << endl << flush;
  }	
}


void Octupole::kick(vektor& R1, vektor& R0, double ds){
  R1[0] = R0[0];
  R1[1] = R0[1]+1.0/6.0*strength_h*ds*(pow(R0[0], 3) - 0.0*3.0*R0[0]*pow(R0[2], 2));
  R1[2] = R0[2];
  R1[3] = R0[3]-1.0/6.0*strength_h*ds*(0.0*3.0*R0[2]*pow(R0[0], 2) - pow(R0[2], 3));
  R1[4] = R0[4];
  R1[5] = R0[5];
}


Chrom::Chrom(double tunex0, double tuney0, double R){
  chrom_x = pow(tunex0/R, 2);
  chrom_y = pow(tuney0/R, 2);
}


void Chrom::kick(vektor& R1, vektor& R0, double ds){
  R1[0] = R0[0];
  R1[1] = R0[1]-chrom_x*ds*R0[5]*R0[0];
  R1[2] = R0[2];
  R1[3] = R0[3]-chrom_y*ds*R0[5]*R0[2];
  R1[4] = R0[4];
  R1[5] = R0[5];	
}


TuneShift::TuneShift(double tunex0, double tuney0, double dtunex,
		     double dtuney, double R){
  coeffx = 2.0*tunex0*dtunex/pow(R, 2);
  coeffy = 2.0*tuney0*dtuney/pow(R, 2);
}


void TuneShift::kick(vektor& R1, vektor& R0, double ds){
  R1[0] = R0[0];
  R1[1] = R0[1]+coeffx*ds*R0[0];
  R1[2] = R0[2];
  R1[3] = R0[3]+coeffy*ds*R0[2];
  R1[4] = R0[4];
  R1[5] = R0[5];	
}


AmplitudeDetuning::AmplitudeDetuning(double tunex0, double tuney0, double
				     Qxs, double Qys, double R, SectorMap& M){
  SMap = M;
  coeffx = tunex0*Qxs/pow(R, 2);  // 2.0 removed
  coeffy = tuney0*Qys/pow(R, 2);  // 2.0 removed	
}


void AmplitudeDetuning::kick(vektor& R1, vektor& R0, double ds){
  emitx = 1.0/SMap.get_betx()*(pow(R0[0], 2)+pow(SMap.get_betx()*R0[1], 2));
  emity = 1.0/SMap.get_bety()*(pow(R0[2], 2)+pow(SMap.get_bety()*R0[3], 2));

  R1[0] = R0[0];
  R1[1] = R0[1]+coeffx*ds*emitx*R0[0];
  R1[2] = R0[2];
  R1[3] = R0[3]+coeffy*ds*emity*R0[2];
  R1[4] = R0[4];
  R1[5] = R0[5];	
}
