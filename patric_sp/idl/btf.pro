pro btfplot, data_dir

!P.MULTI=[0,1,1]
!P.CHARSIZE=1.7
!P.CHARTHICK=3.0
!P.THICK=5.0
!Y.THICK=5.0
!X.THICK=5.0

;Eingabeparameter:

druck=1     ; Ausdruck ja/nein (1/0)

emittance=0
fft_amp=0
fft_phase=0
btf_inv=0

smoothing=20

;IDL device commands:

loadct,12

set_plot,"x"
if(druck EQ 1) then begin
 set_plot,"ps"
 device,filename="btf.eps",/color,/encapsulated,/inches,ysize=4.8;9.6
endif 

; read parameter file idl.dat:

openr, 1, data_dir+"idl.dat"
readf,1,numprocs ; Prozesse
readf,1,e_kin    ; kinetische Energie/MeV
readf,1,qb       ; Ladung/qe
readf,1,mb       ; Masse/mp
readf,1,current  ; Strom      
readf,1,C        ; Ringumfang
readf,1,Nelements    ; elements in cell 
readf,1,cell_length   ; Zellenlaenge
readf,1,Np       ; Teilchen fuer Plot
readf,1,pipe_radius
readf,1,NX
readf,1,NY
close,1

; number of cells in ring:

Ncell=C/cell_length

; coherent (bare) machine tune

baretune=3.24333 ;5.0-Ncell*97.3/360.0

common share2, tunevec2, btf_vec2, btf_amp2, btf_phase2, maxindex2

; Constants:

qe=1.6022e-19
mp=1.6605e-27
I0=current
pi=3.14159265359
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   

; Oeffnen der Datenfiles:

num_p=long(0)
openr, 1, data_dir+"ppic.dat"
tem=fltarr(18)   
while(not EOF(1)) do begin
 readf,1,tem
 num_p=num_p+1
end
close,1

openr, 1, data_dir+"ppic.dat"
ppic=fltarr(18,num_p) 
readf,1,ppic
close,1

; Grid:

s=fltarr(num_p)
s(0:num_p-1)=ppic(0,0:num_p-1)

; Emittance:

ex=fltarr(num_p)
ey=fltarr(num_p)

ex(0:num_p-1)=ppic(1,0:num_p-1)
ey(0:num_p-1)=ppic(2,0:num_p-1)

; offsets

offsetx=fltarr(num_p)
offsety=fltarr(num_p)
pickupx=fltarr(num_p)
pickupy=fltarr(num_p)

offsetx(0:num_p-1)=ppic(13,0:num_p-1)
offsety(0:num_p-1)=ppic(14,0:num_p-1)
pickupx(0:num_p-1)=ppic(16,0:num_p-1)
pickupy(0:num_p-1)=ppic(17,0:num_p-1)


; FFT of offsets:

han=hanning(num_p,/double)

offsetx_fft=fft(han*offsetx,/double)
offsety_fft=fft(han*offsety,/double)
pickupx_fft=fft(han*pickupx,/double)
pickupy_fft=fft(han*pickupy,/double)

Dt=(s(1)-s(0))/(beta0*clight)  ;  time step in track
freq=findgen(num_p/2+1) / (num_p*Dt)  ; resulting frequency grid
tunevec=freq/(beta0*clight/C)

; dipole excitation signal:

dtheta=fltarr(num_p)
dtheta(0:num_p-1)=ppic(15,0:num_p-1)
dtheta_fft=fft(han*dtheta,/double)

; calculate btf from signal 

btf_vec=offsetx_fft/dtheta_fft
btf_amp=abs(btf_vec)
btf_phase=180.0*atan(-dtheta_fft(0:num_p/2)/offsetx_fft(0:num_p/2),/phase)/!PI+90.0
for j=0L, num_p/2 do begin
  if btf_phase[j] GT 180.0 then btf_phase[j]=btf_phase[j]-360.0
endfor

;num_p2=num_p
;tunevec2=tunevec
;btf_amp2=btf_amp
;btf_phase2=btf_phase

; analytic btf (from file):

Nbtf=200
openr,1,'oct.dat'
uv_data=fltarr(3,Nbtf) 
readf,1,uv_data
close,1

btf_analytic=complexarr(Nbtf)
zfreq=fltarr(Nbtf)
Seff=0.012 ; 0.0085

Nbtf=200
openr,1,'octsc.dat'
uv_data_2=fltarr(3,Nbtf) 
readf,1,uv_data_2
close,1

btf_analytic_2=complexarr(Nbtf)

for j=0L, Nbtf-1 do begin
  zfreq[j]=uv_data[0,j]
  btf_analytic[j]=complex(uv_data[1,j],uv_data[2,j])
endfor

for j=0L, Nbtf-1 do begin
  btf_analytic_2[j]=complex(uv_data_2[1,j],uv_data_2[2,j])
endfor

if( fft_amp eq 1) then begin

maxbtf=max(smooth(abs(offsetx_fft),smoothing))
maxindex=!C
maxtune=baretune-0.15  ; tunevec[maxindex] !!!!!!!!!!
;maxbtf=btf_amp[!C]
minbtf=min(smooth(btf_amp,smoothing))
maxana=max(abs(btf_analytic))

print, maxtune

plot,tunevec,smooth(btf_amp,smoothing)^2,/ylog,$
/xstyle,/ystyle,title="",xtitle="tune",ytitle="P [a.u.]",$xrange=[0,3]
xrange=[0.95*maxtune,1.05*maxtune],yrange=[1.0e-3*maxbtf^2,4.0*maxbtf^2]

oplot,tunevec2,smooth(btf_amp2,smoothing)^2,color=190

;plotS, [maxtune,maxtune],[minoffset,maxoffset],line=2
plotS, [baretune-0.15,baretune-0.15],[1.0e-3*maxbtf^2,4.0*maxbtf^2],line=1

;oplot,maxtune+zfreq*Seff,1.0*(maxbtf/maxana)^2*abs(btf_analytic)^2,color=150,line=1
;oplot,maxtune+zfreq*Seff,1.0*(maxbtf/maxana)^2*abs(btf_analytic_2)^2,color=100,line=2

endif

if( fft_phase eq 1) then begin

plot,tunevec,smooth(btf_phase,smoothing),$
/xstyle,/ystyle,title="",ytitle="Phase [deg]",xtitle="tune",xrange=[0.95*maxtune,1.05*maxtune],yrange=[-180,270]

oplot,tunevec2,smooth(btf_phase2,smoothing),color=190

;oplot,maxtune+zfreq*Seff,0.0*90.0+180.0*atan(btf_analytic,/phase)/!PI,color=150,line=1
;oplot,maxtune+zfreq*Seff,0.0*90.0+180.0*atan(btf_analytic_2,/phase)/!PI,color=100,line=2

plotS, [baretune-0.15,baretune-0.15],[-180,270],line=1

endif

;plot emittance:

if( emittance eq 1) then begin

plot,s/C,ex/ex(0),/xstyle,/ystyle,title="",$
xtitle="t [ms]",ytitle="!7e!3!Lx,y!N/!7e!3!L0!N"; /nodata
;oplot,s/C,ex/ex(0),color=150
oplot,s/C,ey/ey(0),color=50

endif

if ( btf_inv eq 1) then begin

maxbtf=max(smooth(abs(offsetx_fft),smoothing))
maxindex=!C
maxtune=baretune  ; tunevec[maxindex] !!!!!!!!!!


btf_vec_real_smooth=smooth(real_part(btf_vec),smoothing)
btf_vec_im_smooth=smooth(imaginary(btf_vec),smoothing)

for j=0L, num_p/2 do begin
  btf_vec[j]=complex(btf_vec_real_smooth[j],btf_vec_im_smooth[j])  
endfor

;maxindex2=maxindex
;btf_vec2=btf_vec

plot,1.0*C/(2.0*!PI*baretune*Seff*beta0*clight)*1.0/smooth(abs(btf_vec[maxindex-100:maxindex+100]),smoothing),$
!PI*smooth(-btf_phase[maxindex-100:maxindex+100],smoothing)/180.0,$
/polar,xrange=[-0.5,1.5],yrange=[-20.0,15.0],/ystyle,ytitle="U",xtitle="V",title=""

oplot,1.0*C/(2.0*!PI*baretune*Seff*beta0*clight)*1.0/smooth(abs(btf_vec2[maxindex2-100:maxindex2+100]),smoothing),$
!PI*smooth(-btf_phase2[maxindex2-100:maxindex2+100],smoothing)/180.0,/polar,color=190

plotS, [0.0,0.0],[-20,15.0],line=1

;oplot,abs(1.0/btf_analytic),(-0.0*0.5*!PI+atan(1.0/btf_analytic,/phase)),/polar,color=150,line=2
;oplot,abs(1.0/btf_analytic_2),-0.0*0.5*!PI+atan(1.0/btf_analytic_2,/phase),/polar,color=100,line=2

endif

if(druck EQ 1) then begin
 device,/close
 spawn,'display btf.eps' 
 set_plot,"x"
endif

END





