#include "obflib.h"
#include "SectorMap.h"
#include "Pic.h"
#include "mathtools.h"

#include <mpi.h>
#include <stdlib.h>


  //-------------------------- Initial distributions ---------------------------

  // -------------------- longitudinal --------------------

void Pic::parabolic_dc(double bunchfactor, double length, double dp0, long Np, long *d,
		       int dummyi){
  length *= 2.0/3.0*bunchfactor;
  dp0 *= sqrt(5.0);
  long i;
  double r1 ,r2;
  coasting_beam(length, Np, d);
  for(i=0; i<pics.size(); ++i){
    r1 = ran1(d);
    r2 = ran1(d);
    pics[i].dp = dp0*sin(2.0*PI*r2) * sqrt(1.0-pow(1.0-r1, 2.0/3.0));
  }
}


void Pic::parabolic(double zlm, double z0, double dp0, long Np, long *d, int dummy){
  dp0 *= sqrt(5.0);
  long i;
  double r1 ,r2, parab_b = 0.0;
  Particle tem;
  for(i=0; i<Np; ++i){
    r1 = ran1(d);
    r2 = ran1(d);
    tem.z = z0 + zlm * sqrt(1.0-pow(1.0-r1, 2.0/3.0))*cos(2.0*PI*r2);
    tem.dp = (dp0*sin(2.0*PI*r2)
	      + parab_b*zlm*cos(2.0*PI*r2))*sqrt(1.0-pow(1.0-r1, 2.0/3.0));
    if(tem.z >= z1 && tem.z < z2)
      pics.push_back(tem);
  }
}


void Pic::coast_gauss(double bunchfactor, double length, double dp0, long Np, long *d,  // SP
		      int dummyi){
  coasting_beam(bunchfactor*length, Np, d);
  gaussz(dp0, d);  // ATTENTION: different 'd' than in original main.cpp
}


void Pic::coasting_beam(double length, long Np, long *d){
  long j;
  Particle tem;
  for(j=0; j<Np; ++j){
    tem.z = length*(ran1(d)-0.5);
    if(tem.z >= z1 && tem.z < z2)
      pics.push_back(tem);
  }
}


void Pic::gaussz(double dp0, long *d){
  long j;
  double u, v, s;
  for(j=0; j<pics.size(); j+= 2){
    do{
      u = 2.0 * ran1(d) - 1.0;
      v = 2.0 * ran1(d) - 1.0;
      s = u*u + v*v;
    }while (s >= 1.0 );

    pics[j].dp   = dp0 * u * sqrt(-2.0*log(s)/s);
    pics[j+1].dp = dp0 * v * sqrt(-2.0*log(s)/s);
  }
}


// Gaussian bunch with cut distribution

void Pic::bunch_gauss(double zlm, double circum, double dp0, long Np, long *d, int dummy){
  zlm /= sqrt(5.0);
  long j;
  int linrf = 1;
  double cut = 8.0;  // 2*sigmas
  Particle tem1, tem2;
  double R = circum/(2.0*PI);
  double Y2u, Y2v, Ym = 1.0-cos(cut*zlm/R), dpu, dpv;
  if(linrf == 1)
    Ym = 0.5*pow(cut*zlm/R, 2);
  double u, v, s, up, vp, sp, zu, zv;
  for(j=0; j<Np; j+= 2){
    do{
      u = 2.0 * ran1(d) - 1.0;
      v = 2.0 * ran1(d) - 1.0;
      s = u*u + v*v;
      up = 2.0 * ran1(d) - 1.0;
      vp = 2.0 * ran1(d) - 1.0;
      sp = up*up + vp*vp;
      zu = zlm * u * sqrt(-2.0*log(s)/s);
      zv = zlm * v * sqrt(-2.0*log(s)/s);
      Y2u = 1.0-cos(zu/R);
      Y2v = 1.0-cos(zv/R);
      dpu = up*sqrt(-2.0*log(sp)/sp);
      dpv = vp*sqrt(-2.0*log(sp)/sp);
      if (linrf == 1){
	Y2u = 0.5*pow(zu/R, 2);
	Y2v = 0.5*pow(zv/R, 2);
      }
    }while(s >= 1.0 || sp >= 1.0 || pow(dpu/cut, 2)+Y2u/Ym > 1.0 ||
	   pow(dpv/cut, 2)+Y2v/Ym > 1.0  );
    tem1.z = zu;
    tem2.z = zv;
    tem1.dp = dpu*dp0;
    tem2.dp = dpv*dp0;
    if(z1 == -zlm*sqrt(5.0) && tem1.z < z2)
      pics.push_back(tem1);
    else if(z2 >= zlm*sqrt(5.0) && tem1.z >= z1)
      pics.push_back(tem1);
    else if(tem1.z >= z1 && tem1.z < z2)
      pics.push_back(tem1);
    if(z1 == -zlm*sqrt(5.0) && tem2.z < z2)
      pics.push_back(tem2);
    else if(z2 >= zlm*sqrt(5.0) && tem2.z >= z1)
      pics.push_back(tem2);
    else if(tem2.z >= z1 && tem2.z < z2)
      pics.push_back(tem2);
  }
}


void Pic::bunch_const(double zlm, double circum, double dp0, long Np, long *d,
		      int linrf){
  dp0 *= sqrt(5.0);
  long i;
  double r1, r2;
  double R = circum/(2.0*PI);
  double Y2, Ym = 1.0-cos(zlm/R);
  if(linrf == 1)
    Ym = 0.5*pow(zlm/R, 2);
  Particle tem;
  for(i = 0; i<Np; ++i){
    do{
      r1 = dp0*(2.0*ran1(d)-1.0);
      r2 = zlm*(2.0*ran1(d)-1.0);
      Y2 = 1.0-cos(r2/R);
      if (linrf == 1) Y2 = 0.5*pow(r2/R, 2);
    }while(pow(r1/dp0, 2)+Y2/Ym > 1.0);
    tem.z = r2;
    tem.dp = r1;
    if(tem.z >= z1 && tem.z < z2)
      pics.push_back(tem);
  }
}


void Pic::barrier_air_bag(double zlm, double dummyd, double dp0, long Np, long *d,
			  int dummyi){
  dp0 *= sqrt(5.0);
  long i;
  double r1, r2;
  Particle tem;
  for(i=0; i<Np; ++i){
    r1 = (2.0*ran1(d)-1.0);
    tem.dp = r1/abs(r1)*dp0;
    tem.z = zlm*(2.0*ran1(d)-1.0);
    if(tem.z >= z1 && tem.z < z2)
      pics.push_back(tem);  	
  }
}


void Pic::bunch_air_bag(double zlm, double circum, double dp0, long Np, long *d,
			int dummyi){
  dp0 *= sqrt(5.0);
  long i;
  double r1, r2;
  double R = circum/(2.0*PI);
  double Y2, rho1, Ym = 1.0-cos(zlm/R);
  Particle tem;
  for(i = 0; i<Np; ++i){
    do{
      r1 = 50.0*(2.0*ran1(d)-1.0);
      r2 = zlm*(2.0*ran1(d)-1.0);
      Y2 = 1.0-cos(r2/R);
      rho1 = 1.0/sqrt(1.0-Y2/Ym);
    }while(rho1<abs(r1));
    tem.z = r2;
    tem.dp = r1/abs(r1)*dp0*sqrt(1.0-Y2/Ym);
    if(tem.z >= z1 && tem.z < z2)
      pics.push_back(tem);
  }
}


  // -------------------- transverse --------------------

void Pic::waterbag_xy(double emittance_x, double emittance_y, double alpha_x,
		      double alpha_y, double beta_x, double beta_y, double D0,
		      double Ds0, double centering, double inj_angle, long *d){
  emittance_x *= 6.;
  emittance_y *= 6.;
  long j;
  double x, y, xs, ys;
  double xmax = sqrt(emittance_x*beta_x);
  double xsmax = sqrt(emittance_x * (1.0+alpha_x*alpha_x) / beta_x);
  double ymax = sqrt(emittance_y*beta_y);
  double ysmax = sqrt(emittance_y * (1.0+alpha_y*alpha_y) / beta_y);

  for(j=0; j<pics.size(); ++j){
    do{
      x = xmax*(2.0*ran1(d)-1.0);
      y = ymax*(2.0*ran1(d)-1.0);
      xs = xsmax*(2.0*ran1(d)-1.0);
      ys = ysmax*(2.0*ran1(d)-1.0);
  }while((x*x/beta_x + beta_x * pow(xs+alpha_x/beta_x*x, 2)) / emittance_x
	   + (y*y/beta_y + beta_y * pow(ys+alpha_y/beta_y*y, 2)) / emittance_y > 1.0);
    pics[j].x = x + D0*pics[j].dp + centering;  // stuff for head tail modes removed; SP
    pics[j].xs = xs + Ds0*pics[j].dp + inj_angle;  // stuff for head tail modes removed; SP
    pics[j].y = y;
    pics[j].ys = ys;
  }
  init_old_coord();
}


void Pic::KV_xy(double emittance_x, double emittance_y, double alpha_x, double alpha_y,
		double beta_x, double beta_y, double D0, double Ds0, double centering,
		double inj_angle, long *d){
  emittance_x *= 4.;
  emittance_y *= 4.;
  long j;
  double u, v, w, rx, rsx, ry, rsy;

  for(j=0; j<pics.size(); ++j){
    u = ran1(d);
    v = ran1(d);
    w = ran1(d);

    rx = sqrt(u*emittance_x*beta_x);
    ry = sqrt((1.-u)*emittance_y*beta_y);
    rsx = rx/beta_x;
    rsy = ry/beta_y;

    pics[j].x = rx * cos(2.*PI*v) + D0*pics[j].dp + centering;
    pics[j].y = ry * cos(2.*PI*w);
    pics[j].xs = -rsx * (alpha_x*cos(2.*PI*v) + sin(2.*PI*v)) + Ds0*pics[j].dp + inj_angle;
    pics[j].ys = -rsy * (alpha_y*cos(2.*PI*w) + sin(2.*PI*w));
  }
  init_old_coord();
}


void Pic::SG(double emittance_x, double emittance_y, double alpha_x,
	     double alpha_y, double beta_x, double beta_y, double D0,
	     double Ds0, double centering, double inj_angle, long *d){
  long j;
  double u, v, w, s, rx, ry;

  double xsrms = sqrt(emittance_x*(1.0+pow(alpha_x, 2))/beta_x);
  double ysrms = sqrt(emittance_y*(1.0+pow(alpha_y, 2))/beta_y);
  double xmax = 4.*emittance_x*beta_x;
  double ymax = 4.*emittance_y*beta_y;

  for(j=0; j<pics.size(); ++j){
    u = ran1(d);
    v = ran1(d);
    w = ran1(d);

    rx = sqrt(u*xmax);
    ry = sqrt(u*ymax);
    pics[j].x = rx * cos(2.0*PI*v) + D0*pics[j].dp + centering;
    pics[j].y = ry * sin(2.0*PI*v);

    do{
      u = 2.0 * ran1(d) - 1.0;
      v = 2.0 * ran1(d) - 1.0;
      s = u*u + v*v;
    }while (s >= 1.0);

    pics[j].xs = xsrms * u * sqrt(-2.0*log(s)/s) + Ds0*pics[j].dp + inj_angle;
    pics[j].ys = ysrms * v * sqrt(-2.0*log(s)/s);
  }
  init_old_coord();
}


void Pic::Gauss_xy(double emittance_x, double emittance_y, double alpha_x,
		   double alpha_y, double beta_x, double beta_y, double D0,
		   double Ds0, double centering, double inj_angle, long *d){
  long j;
  double u, v, w, s, rx, ry, rxs, rys;

  for(j=0; j<pics.size(); ++j){
    do{
      u = 2.0 * ran1(d) - 1.0;
      v = 2.0 * ran1(d) - 1.0;
      w = u*u + v*v;
    }while (w >= 1.0);
    rx = sqrt(emittance_x*beta_x) * u * sqrt(-2.0*log(w)/w);
    ry = sqrt(emittance_y*beta_y) * v * sqrt(-2.0*log(w)/w);

    do{
      u = 2.0 * ran1(d) - 1.0;
      v = 2.0 * ran1(d) - 1.0;
      w = u*u + v*v;
    }while (w >= 1.0);
    rxs = sqrt(emittance_x/beta_x) * u * sqrt(-2.0*log(w)/w);
    rys = sqrt(emittance_y/beta_y) * v * sqrt(-2.0*log(w)/w);

    pics[j].x = rx + D0*pics[j].dp + centering;
    pics[j].y = ry;
    pics[j].xs = rxs - alpha_x/beta_x*rx + Ds0*pics[j].dp + inj_angle;
    pics[j].ys = rys - alpha_y/beta_y*ry;
  }
  init_old_coord();
}


// ---------------------- Teilchenbewegung ---------------------------

void Pic::transport(SectorMap* M, double K_xix, double K_xiy, double boundary){
  vektor R0(6), R1(6);

  for(long j=0; j<pics.size(); ++j){
    if(j >= pics.size())
      break;  // j not always < pics.size() ?; SP
    R0[0] = pics[j].x;
    R0[1] = pics[j].xs;
    R0[2] = pics[j].y;
    R0[3] = pics[j].ys;
    R0[4] = pics[j].z/SP->beta0;
    R0[5] = SP->beta0*pics[j].dp;
    M->transport(R1, R0, K_xix, K_xiy);  // why function call instead of in place evaluation (->efficiency)?; SP
    pics[j].x = R1[0];
    pics[j].xs = R1[1];
    pics[j].y = R1[2];
    pics[j].ys = R1[3];
  //pics[j].z = R1[4]*SP->beta0;
    pics[j].z -= SP->eta0*pics[j].dp*M->get_L();
  //pics[j].dp = R1[5]/SP->beta0;

    if(fabs(pics[j].x) >= boundary || fabs(pics[j].y) >= boundary)
      pics.erase(pics.begin()+j);
  }
}

long Pic::localLoss_x(double loBound, double upBound){
  long lolo = 0;
  for(long j=0; j<pics.size(); ++j){
    if(pics[j].x <= loBound || pics[j].x >= upBound){
      pics.erase(pics.begin()+j);
      ++lolo;
    }
  }
  return lolo;
}


void Pic::kick(Grid2D& Ex, Grid2D& Ey, double ds){
  double beta0 = SP->beta0, gamma0 = SP->gamma0, A = SP->A, Z = SP->Z;
  double temp = qe*Z/(mp*A*pow(gamma0, 3)*pow(beta0*clight, 2))*ds;
  double ex, ey;

  for(long j=0; j < pics.size(); ++j){
    ex = temp*Ex.Grid2PIC(pics[j].x, pics[j].y);
    ey = temp*Ey.Grid2PIC(pics[j].x, pics[j].y);

  //momenta:

    pics[j].xs += ex;
    pics[j].ys += ey;
  }
}


void Pic::kick(Grid3D& Ex, Grid3D& Ey, double ds){
  double beta0 = SP->beta0, gamma0 = SP->gamma0, A = SP->A, Z = SP->Z;
  double temp = qe*Z/(mp*A*pow(gamma0, 3)*pow(beta0*clight, 2))*ds;
  double ex, ey;

  for(long j=0; j<pics.size(); ++j){
    ex = temp*Ex.Grid2PIC(pics[j].x, pics[j].y, pics[j].z);
    ey = temp*Ey.Grid2PIC(pics[j].x, pics[j].y, pics[j].z);

  //momenta:

    pics[j].xs += ex;
    pics[j].ys += ey;
  }
}


void Pic::kick(double fx, double fy){
  for(long j=0; j<pics.size(); ++j){
    pics[j].xs += fx;
    pics[j].ys += fy;
  }
}


void Pic::cavity_kick(double voltage0, int harmonic, double R){
  double beta0 = SP->beta0, gamma0 = SP->gamma0, A = SP->A, Z = SP->Z;
  for(long j=0; j<pics.size(); ++j)
    pics[j].dp -= Z*qe*voltage0*sin(harmonic*pics[j].z/R)/(A*mp*gamma0*pow(beta0*clight, 2));
}


void Pic::cavity_kick_linear(double voltage0, int harmonic, double R){
  double beta0 = SP->beta0, gamma0 = SP->gamma0, A = SP->A, Z = SP->Z;
  for(long j=0; j<pics.size(); ++j)
    pics[j].dp -= Z*qe*voltage0*(harmonic*pics[j].z/R)/(A*mp*gamma0*pow(beta0*clight, 2));
}


void Pic::barrier_kick(double zm1, double zm2){
  for(long j=0; j<pics.size(); ++j){
    if (pics[j].z <= zm1 && pics[j].dp < 0.0)
      pics[j].dp =-pics[j].dp;
    if (pics[j].z >= zm2 && pics[j].dp > 0.0)
      pics[j].dp =-pics[j].dp;	
  }
}


void Pic::kick(ThinLens& M, double ds){
  vektor R0(6), R1(6);
  for(long j=0; j<pics.size(); ++j){
    if(j >= pics.size())
      break;
    R0[0] = pics[j].x;
    R0[1] = pics[j].xs;
    R0[2] = pics[j].y;
    R0[3] = pics[j].ys;
    R0[4] = pics[j].z;
    R0[5] = pics[j].dp;
    M.kick(R1, R0, ds);
    pics[j].x = R1[0];
    pics[j].xs = R1[1];
    pics[j].y = R1[2];
    pics[j].ys = R1[3];
    pics[j].z = R1[4];
    pics[j].dp = R1[5];
  }
}


void Pic::impedance_kick(Grid1D& kick, double circum, double ds){
  for(long j=0; j<pics.size(); ++j)
    pics[j].xs += ds/circum*kick.Field2Pic(pics[j].z);
}


void Pic::linear_SC_kick(double dQxm, double dQym, double tunex, double tuney,
                         Grid1D& ldy, double ldy0,
			 Grid1D& dipole_current_x, Grid1D& dipole_current_y,
			 double circum, double ds){
  double R = circum/(2.0*PI);
  double scfact;
  double offsetx;
  double offsety;

  for(long j=0; j<pics.size(); ++j){
    scfact = ldy.Field2Pic(pics[j].z)/ldy0;
    if(scfact > 0.0){
      offsetx = dipole_current_x.Field2Pic(pics[j].z)/(SP->beta0*clight*ldy0*scfact);
      offsety = dipole_current_y.Field2Pic(pics[j].z)/(SP->beta0*clight*ldy0*scfact);
      pics[j].xs-= scfact*2.0*ds*tunex*dQxm/pow(R, 2)*(pics[j].x-offsetx);
      pics[j].ys-= scfact*2.0*ds*tuney*dQym/pow(R, 2)*(pics[j].y-offsety);
    }
  }
}


void Pic::nonlinear_SC_kick(double xrms, double yrms, double dQxm, double dQym,
			    double tunex, double tuney, Grid1D& ldy, double ldy0,
			    double circum, double ds){
  double R = circum/(2.0*PI);
  double scfact;
  double offsetx = offset_x();
  double offsety = offset_y();

  for(long j=0; j<pics.size(); ++j){
    scfact = ldy.Field2Pic(pics[j].z)/ldy0;

  //if( pow(pics[j].x-offsetx, 2)+pow(pics[j].y-offsety, 2) < 6.0*(pow(xrms, 2)+pow(yrms, 2)) )
  //{
    pics[j].xs-= scfact*2.0*ds*tunex*dQxm/pow(R, 2)*
      (pics[j].x-offsetx)*
      (1.0-1.0/18.0*(2.0*xrms+yrms)/(xrms*xrms*(xrms+yrms))*pow(pics[j].x-offsetx, 2)-1.0/(6.0*yrms*(xrms+yrms))*pow(pics[j].y-offsety, 2));
    pics[j].ys-= scfact*2.0*ds*tuney*dQym/pow(R, 2)*
      (pics[j].y-offsety)*(1.0-1.0/18.0*(2.0*yrms+xrms)/(yrms*yrms*(xrms+yrms))
			   *pow(pics[j].y-offsety, 2)-1.0/(6.0*xrms*(xrms+yrms))*pow(pics[j].x-offsetx, 2));
  //}
  }
}


void Pic::dipole_kick_simple(double dQxm, double dQym, double tunex, double tuney,
			     Grid1D& ldy, double ldy0, Grid1D& dipole_current_x,
			     Grid1D& dipole_current_y, double circum, double ds){
  double R = circum/(2.0*PI);
  double scfact;
  double offsetx;
  double offsety;

  for(long j=0; j<pics.size(); ++j){
    scfact = ldy.Field2Pic(pics[j].z)/ldy0;
    if( scfact > 0.0 ){
      offsetx = dipole_current_x.Field2Pic(pics[j].z)/(SP->beta0*clight*ldy0*scfact);
      offsety = dipole_current_y.Field2Pic(pics[j].z)/(SP->beta0*clight*ldy0*scfact);
      pics[j].xs-= scfact*2.0*ds*tunex*dQxm/pow(R, 2)*offsetx;
      pics[j].ys-= scfact*2.0*ds*tuney*dQym/pow(R, 2)*offsety;
    }
  }
}


/*!
  Bandwidth limited noise dipole modulation kick for BTF.

  Input

  time: t
  step: ds
  kick angle: theta
  bandwidth: freq0, freq1
  simulation duration: tend
  harmonic: n
  random number initialization: d

  Output

  xs kick
*/


double Pic::dipole_mod_kick(double t, double ds, double circum, double theta, double freq0,
	                    double freq1, double tend, int n, long* d){
  double dtheta = 0.0, R = circum/(2.0*PI), beta0 = SP->beta0;
  double dfreq = 1.0/(tend);
  const int Nran_max = 10000;
  int Nran = (int)floor((freq1-freq0)/dfreq);
  static int flag = 0;
  static double freq[Nran_max];
  static double phase[Nran_max];
  if(Nran >= Nran_max){
    cout << "Dipole kick: Nran > Nran_max!\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  if(flag == 0){
    for(int j=0; j<Nran; ++j){
  //freq[j] = freq0+(freq1-freq0)*ran1(d);
      freq[j] = floor(freq0/dfreq)*dfreq+j*dfreq;
      phase[j] = PI*ran1(d);
    }	
    flag = 1;	
  }

  for(int l=0; l<Nran; ++l)
    dtheta += theta*cos(2.0*PI*freq[l]*t+phase[l])/Nran;

  if(t<20.0/freq0)
    dtheta = dtheta*t*freq0/20.0;

  for(long j=0; j<pics.size(); ++j){
    pics[j].xs += ds/R*dtheta*cos(n/R*pics[j].z);
  }

  return 0.5*dtheta*beta0*clight*cos(n/R*beta0*clight*t);
}


double Pic::dipole_mod_kick(double t, double ds, double circum, double theta,
			    double freq, int n){
  double dtheta = 0.0, R = circum/(2.0*PI), beta0 = SP->beta0;

  dtheta = theta*cos(2.0*PI*freq*t);

  if(t<20.0/freq)
    dtheta = dtheta*t*freq/20.0;

  for(long j=0; j<pics.size(); ++j){
    pics[j].xs += ds/R*dtheta*sin(n/R*pics[j].z);
  }

  return 0.5*dtheta*beta0*clight*sin(n/R*beta0*clight*t);
}



double Pic::pickup_signal(Grid1D& dipole_current, double circum, double t){
  double beta0 = SP->beta0;
  double z0 = 0.5*circum-fmod(clight*beta0*t, circum);  // -0.5*circum;
  return dipole_current.get_grid_lin(z0);
}


  //----------------------intrabeam-scattering------------------------------


void Pic::langevin(double beta_fxy, double beta_fz, double Dxy, double Dz,
		   double ds, double beta_x0, double beta_y0, long* d){
  double Rx, Ry, Rz;
  double beta0 = SP->beta0, gamma0 = SP->gamma0, A = SP->A, Z = SP->Z;
  for(long j=0; j<pics.size(); ++j){
    Rx = 0.0; Ry = 0.0; Rz = 0.0;
    for(int l=1; l <= 10; ++l){
      Rx += sqrt(24.0/10.0)*(ran1(d)-0.5);
      Ry += sqrt(24.0/10.0)*(ran1(d)-0.5);
      Rz += sqrt(24.0/10.0)*(ran1(d)-0.5);
    }	

    pics[j].xs += -beta_fxy/(beta0*clight)*pics[j].xs*ds
      +Rx*sqrt(Dxy*ds/(beta_x0*beta0*clight));
    pics[j].x += -beta_fxy/(beta0*clight)*pics[j].x*ds
      +Rx*sqrt(beta_x0*Dxy*ds/(beta0*clight));

    pics[j].ys += -beta_fxy/(beta0*clight)*pics[j].ys*ds
      +Ry*sqrt(Dxy*ds/(beta_y0*beta0*clight));
    pics[j].y += -beta_fxy/(beta0*clight)*pics[j].y*ds
      +Ry*sqrt(beta_y0*Dxy*ds/(beta0*clight));

    pics[j].dp += -beta_fz/(beta0*clight)*pics[j].dp*ds
      +Rz*sqrt(Dz*ds/(beta0*clight));
  }
}

  //--------------------Interpolations: Scatter-Gather---------


void Pic::gatherZ(double pic_charge, Grid1D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Field(pic_charge, pics[j].z);
}


void Pic::gatherX(double pic_charge, Grid1D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Field(pic_charge*pics[j].x, pics[j].z);
}

void Pic::gatherXs(double pic_charge, Grid1D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Field(pic_charge*pics[j].xs, pics[j].z);
}


void Pic::gatherY(double pic_charge, Grid1D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Field(pic_charge*pics[j].y, pics[j].z);
}


void Pic::gatherXY(double pic_charge, Grid2D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Grid(pic_charge, pics[j].x, pics[j].y);
}


void Pic::gatherZX(double pic_charge, Grid2D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Grid(pic_charge, pics[j].z, pics[j].x);
}


void Pic::gatherXXs(double pic_charge, Grid2D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Grid(pic_charge, pics[j].x, pics[j].xs);
}


void Pic::gatherYYs(double pic_charge, Grid2D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Grid(pic_charge, pics[j].y, pics[j].ys);
}


void Pic::gatherXsYs(double pic_charge, Grid2D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Grid(pic_charge, pics[j].xs, pics[j].ys);
}


void Pic::gatherXYZ(double pic_charge, Grid3D& target){
  long j;
  target.reset();
  for(j=0; j<pics.size(); ++j)
    target.Pic2Grid(pic_charge, pics[j].x, pics[j].y, pics[j].z);
}



//--------------------Slice2Slice-----------------------------

vector<Particle> Pic::get_particles_left(double length){
  vector<Particle> tem;
  for(int j=0; j<pics.size(); ++j)
    if(pics[j].z<z1){
      if(pics[j].z < -0.5*length)
	pics[j].z += length;
      tem.push_back(pics[j]);
      pics.erase(pics.begin()+j);
    }	
  return tem;
}


void Pic::periodic_bc(double length){
  for(int j=0; j<pics.size(); ++j){
    if(pics[j].z < -0.5*length)
      pics[j].z += length;
    else
      if(pics[j].z > 0.5*length)
	pics[j].z -= length;
  }
}


vector<Particle> Pic::get_particles_right(double length){
  vector<Particle> tem;
  for(int j=0; j<pics.size(); ++j)
    if(pics[j].z > z2){
      if(pics[j].z > 0.5*length)
	pics[j].z -= length;
      tem.push_back(pics[j]);
      pics.erase(pics.begin()+j);
    }
  return tem;
}


void Pic::add_particles(vector<Particle> &part){
  pics.reserve(pics.size() + part.size());  // more efficient; SP
  pics.insert(pics.end(), part.begin(), part.end());
  //  for(int j=0; j<part.size(); ++j)
  //    pics.push_back(part[j]);
}


void Pic::remove_lost_particles(double length){
  for(int j=0; j<pics.size(); ++j)
    if(pics[j].z < -0.5*length || pics[j].z > 0.5*length)
      pics.erase(pics.begin()+j);
}


  //--------------------Ausgabe----------------------------------


void Pic::print(int subset){
  long j, jran, d = 100, np = pics.size();
  if(subset>np){
    printf("subset > NPIC: Execution aborted.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  float *tem = new float[8];
  for(j=0; j<subset; ++j){
    jran = (int)(ran1(&d)*(np-2.0));  // why -2 instead of -1?; SP
    tem[0] = pics[jran].x;
    tem[1] = pics[jran].xs;
    tem[2] = pics[jran].y;
    tem[3] = pics[jran].ys;
    tem[4] = pics[jran].dp;
    tem[5] = pics[jran].z;
    tem[6] = get_phaseadvance_h(jran);  //   get_wavelength_h(jran);
    tem[7] = get_phaseadvance_v(jran);  //     get_wavelength_v(jran);
    fwrite(tem, sizeof(float), 8, out);
  }
  fflush(out);
  delete tem;
}


  //--------------------Momente-------------------------------------

double Pic::rms_emittance_x(){
  long j, n = pics.size();
  double tem1 = 0.0, tem2 = 0.0, tem3 = 0.0, tem4 = 0.0, tem5 = 0.0;
  for(j=0; j<n; ++j){
    tem4 += pics[j].x;
    tem5 += pics[j].xs;
  }
  tem4 /= n;
  tem5 /= n;

  for(j=0; j<n; ++j){
    tem1 += pow(pics[j].x-tem4, 2);
    tem2 += pow(pics[j].xs-tem5, 2);
    tem3 += (pics[j].x-tem4)*(pics[j].xs-tem5);
  }
  return sqrt(tem1*tem2/pow((double)n, 2)-pow(tem3/n, 2));
}


double Pic::rms_emittance_y(){
  long j, n = pics.size();
  double tem1 = 0.0, tem2 = 0.0, tem3 = 0.0;
  for(j=0; j<n; ++j){
    tem1 += pow(pics[j].y, 2);
    tem2 += pow(pics[j].ys, 2);
    tem3 += pics[j].y*pics[j].ys;
  }
  return sqrt(tem1*tem2/pow((double)n, 2)-pow(tem3/n, 2));
}


double Pic::x_rms(){
  long j, n = pics.size();
  double tem1 = 0.0, tem2 = 0.0;
  for(j=0; j<n; ++j){
    tem1 += pow(pics[j].x, 2);
    tem2 += pics[j].x;
  }
  return sqrt(tem1/pics.size()-pow(tem2/pics.size(), 2));
}


double Pic::offset_x(){
  long j, n = pics.size();
  double tem1 = 0.0;
  for(j=0; j<n; ++j)
    tem1 += pics[j].x;
  return tem1/pics.size();
}


double Pic::offset_y(){
  long j, n = pics.size();
  double tem1 = 0.0;
  for(j=0; j<n; ++j)
    tem1 += pics[j].y;
  return tem1/pics.size();
}


double Pic::y_rms(){
  long j;
  double tem1 = 0.0, tem2 = 0.0;
  for(j=0; j<pics.size(); ++j){
    tem1 += pow(pics[j].y, 2);
    tem2 += pics[j].y;
  }
  return sqrt(tem1/pics.size()-pow(tem2/pics.size(), 2));
}


double Pic::x2y2(){
  long j;
  double tem1 = 0.0, tem2 = 0.0;
  for(j=0; j<pics.size(); ++j){
    tem1 += pow(pics[j].x*pics[j].y, 2);
    tem2 += pics[j].x*pics[j].y;
  }
  return tem1/pics.size()-pow(tem2/pics.size(), 2);
}


double Pic::xy(){
  long j;
  double tem = 0.0;
  for(j=0; j<pics.size(); ++j)
    tem += pics[j].x*pics[j].y;
  return tem/pics.size();
}


double Pic::xzn(double k_mode, double zlm){
  long j;
  double tem = 0.0;
  for(j=0; j<pics.size(); ++j)
    tem += pics[j].x*cos(0.5*k_mode*PI*(pics[j].z+zlm)/zlm);
  return tem/pics.size();
}


double Pic::x_max(){
  long j;
  double tem = fabs(pics[0].x);
  for(j=1; j<pics.size(); ++j){
    if(fabs(pics[j].x) > tem)
      tem = fabs(pics[j].x);
  }
  return tem;
}


double Pic::y_max(){
  long j;
  double tem = fabs(pics[0].y);
  for(j=1; j<pics.size(); ++j){
    if(fabs(pics[j].y) > tem)
      tem = fabs(pics[j].y);
  }
  return tem;
}


double Pic::z_max(){
  long j;
  double tem = pics[0].z;
  for(j=1; j<pics.size(); ++j){
    if(pics[j].z > tem)
      tem = pics[j].z;
  }
  return tem;
}


double Pic::z_min(){
  long j;
  double tem = pics[0].z;
  for(j = 1; j<pics.size(); ++j){
    if(pics[j].z < tem)
      tem = pics[j].z;
  }
  return tem;
}


double Pic::z_mean(){
  long j;
  double tem = 0.0;
  for(j=0; j<pics.size(); ++j)
    tem += pics[j].z;
  return tem/pics.size();
}


double Pic::z2_mean(){
  long j;
  double tem = 0.0;
  for(j=0; j<pics.size(); ++j)
    tem += pow(pics[j].z, 2);
  return tem/pics.size();
}


double Pic::rms_z_width(){
  return sqrt(fabs(z2_mean()-pow(z_mean(), 2)));
}


double Pic::pz_mean(){
  long j;
  double tem = 0.0;
  for(j=0; j<pics.size(); ++j)
    tem += pics[j].dp;
  return tem/pics.size();
}

double Pic::pz2_mean(){
  long j;
  double tem = 0.0;
  for(j=0; j<pics.size(); ++j)
    tem += pow(pics[j].dp, 2);
  return tem/pics.size();
}


double Pic::rms_momentum_spread(){
  return sqrt(fabs(pz2_mean()-pow(pz_mean(), 2)));
}


double Pic::entropy(Grid2D& target){
  double Hfunc = 0.0;
  int NPIC = get_size();
  int Nxs = target.get_NX();
  int Nys = target.get_NY();
  double dxs = target.get_dx();
  double dys = target.get_dy();
  gatherXsYs(1.0/(double)NPIC, target);
  for(int j=0; j<Nxs; ++j)
    for(int l=0; l<Nys; ++l){
      if(target(j, l) > 0.0)	
	Hfunc += target(j, l)*log(target(j, l))*dxs*dys;
    }
  return Hfunc;
}


  //--------------------Phase advance-------------------------------------

/*
  void Pic::update_wavelength_h(double ds, double offset_x)
  {
  int n = pics.size();
  for(int j=0; j<n; ++j)	
  {
  if(pics[j].x > offset_x && pics[j].x1 > offset_x)
  pics[j].lambda_tmp_h += ds;
  else if(pics[j].x < offset_x && pics[j].x1 < offset_x)
  pics[j].lambda_tmp_h += ds; 	
  else	
  {
  pics[j].lambda_tmp_h += (offset_x-pics[j].x1)*ds/(pics[j].x-pics[j].x1);
  pics[j].lambda_h = pics[j].lambda_tmp_h;
  pics[j].lambda_tmp_h = 0.0;
  }

  pics[j].x1 = pics[j].x;
  }
  }

  void Pic::update_wavelength_v(double ds)
  {
  int n = pics.size();
  for(int j=0; j<n; ++j)	
  {
  if(pics[j].y > 0.0 && pics[j].y1 > 0.0)
  pics[j].lambda_tmp_v += ds;
  else if( pics[j].y < 0.0 && pics[j].y1 < 0.0)
  pics[j].lambda_tmp_v += ds; 	
  else	
  {
  pics[j].lambda_tmp_v += -pics[j].y1*ds/(pics[j].y-pics[j].y1);
  pics[j].lambda_v = pics[j].lambda_tmp_v;
  pics[j].lambda_tmp_v = 0.0;
  }

  pics[j].y1 = pics[j].y;
  }
  }


*/

/*!
  In order to calculate the instantanous
  phase advance we need the exit coordinates
  from the previous two cells
*/


void Pic::init_old_coord(){
  int n = pics.size();
  for(int j=0; j<n; ++j){
    pics[j].x2 = pics[j].xs2 = pics[j].x1 = pics[j].xs1 = 0.;
    pics[j].y2 = pics[j].ys2 = pics[j].y1 = pics[j].ys1 = 0.;
  }
}


void Pic::store_old_coordinates(){
  int n = pics.size();
  for(int j=0; j<n; ++j){
    pics[j].x2 = pics[j].x1;
    pics[j].xs2 = pics[j].xs1;
    pics[j].x1 = pics[j].x;
    pics[j].xs1 = pics[j].xs;
    pics[j].y2 = pics[j].y1;
    pics[j].ys2 = pics[j].ys1;
    pics[j].y1 = pics[j].y;
    pics[j].ys1 = pics[j].ys;
  }
}


/*!
  Calculates the horizontal phase advance per cell from the exit coordinates
  of the previous two cells.
  Does not work with finite dispersion yet.
*/

  //double Pic::get_wavelength_h(int j)
  //{
  // return pics[j].lambda_h;
  //}


double Pic::get_phaseadvance_h(int j){
  double ah, dh;	
  double x = pics[j].x, xs = pics[j].xs, x1 = pics[j].x1, xs1 = pics[j].xs1,
    x2 = pics[j].x2, xs2 = pics[j].xs2;

  if(fabs(x1-x2*xs1/xs2) > 0.0)
    ah = (x-x1*xs1/xs2)/(x1-x2*xs1/xs2);
  else return 0.0;

  if(fabs(xs1-xs2*x1/x2) > 0.0 )
    dh = (xs-xs1*x1/x2)/(xs1-xs2*x1/x2);
  else return 0.0;

  if(fabs(ah+dh)<2.0)
    return acos(0.5*(ah+dh));
  else return 0.0;
}


/*!
  Calculates the vertical phase advance per cell from the exit coordinates
  of the previous two cells.
  Does not work with finite dispersion yet.
*/

  //double Pic::get_wavelength_v(int j)
  //{
  // return pics[j].lambda_v;
  //}


double Pic::get_phaseadvance_v(int j){
  double av, dv;	
  double y = pics[j].y, ys = pics[j].ys, y1 = pics[j].y1, ys1 = pics[j].ys1,
    y2 = pics[j].y2, ys2 = pics[j].ys2;

  if(fabs(y1-y2*ys1/ys2) > 0.0)
    av = (y-y1*ys1/ys2)/(y1-y2*ys1/ys2);
  else return 0.0;

  if(fabs(ys1-ys2*y1/y2) > 0.0 )
    dv = (ys-ys1*y1/y2)/(ys1-ys2*y1/y2);
  else return 0.0;

  if(fabs(av+dv)<2.0)
    return acos(0.5*(av+dv));
  else
    return 0.0;
}

/*
  double Pic::rms_wavelength_h(){
  double tem = 0.0;
  int n = pics.size();
  for(int j=0; j<n; ++j)
  tem += get_wavelength_h(j)/n;
  return tem;
  }


  double Pic::rms_wavelength_v(){
  double tem = 0.0;
  int n = pics.size();
  for(int j=0; j<n; ++j)
  tem += get_wavelength_v(j)/n;
  return tem;
  }

*/


double Pic::rms_phaseadvance_h(){
  double tem = 0.0;
  int n = pics.size();
  for(int j=0; j<n; ++j)
    tem += get_phaseadvance_h(j);
  return tem/double(n);
}


double Pic::rms_phaseadvance_v(){
  double tem = 0.0;
  int n = pics.size();
  for(int j=0; j<n; ++j)
    tem += get_phaseadvance_v(j);
  return tem/double(n);
}
