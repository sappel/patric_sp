#include "obflib.h"
#include "mathtools.h"
#include "TImpedance.h"
#include <algorithm>

/*!
Transverse wake function for a RLC circuit (Ng, p. ).
*/

double Trans_Wake(double z, double omega0, double R, double tune0, double nres, double Rs, double Qs)
{
  double omegar=(nres-tune0)*omega0; 
  double alpha=omegar/(2.0*Qs);
  double omegas=sqrt(pow(omegar,2)-pow(alpha,2));
  double v0=omega0*R;
  double R1s=Rs*omegar/v0;
  double wake=-R1s*v0*omegar/(Qs*omegas)*exp(-alpha*z/v0)*sin(omegas*z/v0);
  if (z > 0) return wake;
  else return 0.0;
}

/*!
Transverse wake function for thick resistive wall (Ng, p. ).
*/


double Trans_Wake(double z, double beta0, double circum, double b, double leit)
{
 double Z0=377.0;	
 double wake=-pow(beta0,1.5)*clight*circum/(PI*pow(b,3)*sqrt(z))*sqrt(Z0/(PI*leit));  
 if (z > 0) return wake;
 else return 0.0;	
}	


/*!
Kick vector in z for fixed time induced by the wake field.
Output: kick vector (should be set to zero !).
*/

void InducedWakeKick(Grid1D& kick,Grid1D& dipole_current,Grid1D& dipole_current_s,
                 double tune, double omega0, double nres, double Rs, double Qs, double b, double leit, 
		 double beta0, double E0, double charge) 
{
  double dz=dipole_current.get_dz();
  int n=dipole_current.get_size();
  double circum=dz*n;
  double R=circum/(2.0*PI);
  double tmp=0.0, psi=0.0, psis=0.0;
  int turns=2; // number of turns considered for integration (max. 3)
  for(int j=0; j<n; j++)    
     {
       tmp=0.0;    
       for(int k=0; k<n*turns; k++)   // integration along z
	    {
	     if (j+k < n)
	            {
		     psi=dipole_current[j+k]/(beta0*clight);  // psi: dipole moment at position z=(j+k)*dz
	             psis=dipole_current_s[j+k]/(beta0*clight);
		    }
	     else if (j+k < 2*n)
		   {
		     psi=dipole_current[j+k-n]/(beta0*clight);  // psi: dipole moment at position z=(j+k)*dz
	             psis=dipole_current_s[j+k-n]/(beta0*clight);
		    }
	     else if (j+k < 3*n)
		   {
		     psi=dipole_current[j+k-2*n]/(beta0*clight);  // psi: dipole moment at position z=(j+k)*dz
	             psis=dipole_current_s[j+k-2*n]/(beta0*clight);
		   }
	     else
		  { psi=0.0; psis=0.0; }
	  
	      tmp+=dz*(psi*cos(-tune*k*dz/R)+R/tune*psis*sin(-tune*k*dz/R))
	           *(Trans_Wake(k*dz,omega0,circum/(2.0*PI),tune,nres,Rs,Qs)
		   + Trans_Wake(k*dz,beta0,circum,b,leit));         
	    }
    
       kick[j]+=-charge/(pow(beta0,2)*E0)*tmp;    
     }
}


/*!
Transverse impedance function.
*/

komplex Trans_Impedance(double omega,double omega0, double tune0, double nres, double Rs, double Qs)
{
   komplex i(0.0,1.0);
   double omegar=(nres-tune0)*omega0; 
   komplex cavimp=(omegar/omega)*Rs/(1.0-i*Qs*(omega/omegar-omegar/omega));
   return cavimp;
}


void InducedKick(Grid1D& kick, Grid1D& dipole_current,   
                        double Zimage, double beta0, double E0, double charge) 
{
 int n=dipole_current.get_size();	
 for(int j=0; j<n; j++)
    kick[j]+=-charge/(beta0*E0)*Zimage*dipole_current[j];       
}


void InducedKick(Grid1D& kick, Grid1D& dipole_current, Grid1D& dipole_current_s,  
                 double tune, double omega0, double nres, double Rs, double Qs, 
		 double beta0, double E0, double charge, double t) 
{
	komplex i(0.0,1.0);
        double bl_amp, bl_real, bl_imag, bl_phase; 
	double dz=dipole_current.get_dz();
	int n=dipole_current.get_size();
	double circum=dz*n;
	vektor psi_1(n), psi_2(n);
	for(int j=0; j<n; j++) 
	{ 
	 psi_1[j]=dipole_current[j]*cos(omega0*tune*t)-circum/(2.0*PI*tune)*dipole_current_s[j]*sin(omega0*tune*t);
	 psi_2[j]=dipole_current[j]*sin(omega0*tune*t)+circum/(2.0*PI*tune)*dipole_current_s[j]*cos(omega0*tune*t);
	}
	
	realft(psi_1,1);
	realft(psi_2,1);
	
	for(int h=0; h < n/2-1; h++)
	{
		komplex temc1(psi_1[2*h],psi_1[2*h+1]);   // dipole current amplitude 
		komplex temc2(psi_2[2*h],psi_2[2*h+1]);   
		if ( h==0) { temc1=komplex(psi_1[0],0.0); 
		             temc2=komplex(psi_2[0],0.0); }
		
		komplex timpp=Trans_Impedance(omega0*(h+tune),omega0,tune,nres,Rs,Qs);// fast wave
		komplex timpm=Trans_Impedance(omega0*(h-tune),omega0,tune,nres,Rs,Qs); // slow wave
		
		komplex temcp=i*timpp*(temc1-i*temc2)*exp(i*tune*omega0*t);
		komplex temcm=i*timpm*(temc1+i*temc2)*exp(-i*tune*omega0*t);
                
		komplex temctot=temcp+temcm;
		
		psi_1[2*h]=temctot.real();
		psi_1[2*h+1]=temctot.imag();

	}

	realft(psi_1,-1);

	for(int j=0; j<n; j++)
		kick[j]+=0.5*charge/(beta0*E0)*psi_1[j];    

	return;
}




void InducedKick(Grid1D& kick,Grid1D& dipole_current,
                 double tune, double omega0, double nres, double Rs, double Qs, 
		 double beta0, double E0, double charge) 
{
	komplex i(0.0,1.0);
        double bl_amp, bl_real, bl_imag, bl_phase; 
	double dz=dipole_current.get_dz();
	int n=dipole_current.get_size();
	double circum=dz*n;
	vektor current(n);
	for(int j=0; j<n; j++) current[j]=dipole_current[j]; 

	realft(current,1);
	
	for(int h=0; h < n/2-1; h++)
	{
		komplex temc(current[2*h],current[2*h+1]);   // dipole current amplitude 
		komplex temcc(current[2*h],-current[2*h+1]); // complex conj.
		if ( h==0) { temc=komplex(current[0],0.0); 
		             temcc=komplex(current[0],0.0); }
		
		komplex timpp=Trans_Impedance(omega0*(h+tune),omega0,tune,nres,Rs,Qs);// fast wave
		komplex timpm=Trans_Impedance(omega0*(h-tune),omega0,tune,nres,Rs,Qs); // slow wave
		
		komplex temc2(0.0,0.0), temc3(0.0,0.0);
		//temc2=i*temc*i*timpm.imag(); // slow wave
		//temc3=i*temc*i*timpp.imag();        // fast wave

		temc2=i*temc*(timpm.real()+i*0.5*(timpm.imag()+timpp.imag())); // slow wave
		temc3=i*temcc*(timpp.real());        // fast wave
		
		bl_real = temc2.real()+temc3.real();
		bl_imag = temc2.imag()+temc3.imag();

		current[2*h]=bl_real;
		current[2*h+1]=bl_imag;

	}

	realft(current,-1);

	for(int j=0; j<n; j++)
		kick[j]+=charge/(beta0*E0)*current[j]; // factor 0.5 added !    

	return;
}



/*! 
Coasting beam (single slice) approximation: 
The induced xprime kick per turn.
offset: the actual <x>(t) in m
current0: the dc current in A
dt: simulation time step ds/(beta0*clight)
E0: total ion energy/u 
charge: ion charge Z*qe
*/

double InducedKick(double offset, double ds, komplex dqc, double beta0, double tune, double circum)
{	
 double dqr=dqc.real(), dqi=dqc.imag();
 static double offset_0=0.0;
 double offset_d=(offset-offset_0)/ds;
 offset_0=offset;
 double omega0=2.0*PI*beta0*clight/circum;
 return 8.0*PI*PI*tune/circum*(-dqr*offset+dqi*beta0*clight*offset_d/(tune*omega0));    
}




