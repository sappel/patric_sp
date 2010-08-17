function Qshift, dQsc, Qs, m
 dQm=0.0	
 if (m gt 0.0) then dQm=-dQsc/2.0+sqrt((dQsc/2.0)^2+(m*Qs)^2) else dQm=-dQsc/2.0-sqrt((dQsc/2.0)^2+(m*Qs)^2)
 return, dQm
end

function fdist, tune, Seff
	fd=0.0
	if (1.0-tune^2/Seff^2 gt 0.0) then fd=(1.0-tune^2/Seff^2)^2 
	return, fd
end	

pro schottky, data_dir

!Path='/u/boine/codes/idllib:' + !Path
!Path='/u/boine/codes/idllib/textoidl:' + !Path

!P.MULTI=[0,1,1]
!P.CHARSIZE=1.7
!P.CHARTHICK=3.0
!P.THICK=5.0
!Y.THICK=5.0
!X.THICK=5.0

;Eingabeparameter:

druck=1    ; Ausdruck ja/nein (1/0)
nmax=0     ; 2*4096+2*1024   ; number of cells used for fft
nwin=1     ; Windows
offset=0
fft_offset=1


;IDL device commands:

loadct,12

set_plot,"x"
if(druck EQ 1) then begin
 set_plot,"ps"
 device,filename="schottky.eps",/color,/encapsulated,/inches,ysize=4.8 ; 9.6
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

; Constants:

qe=1.6022e-19
mp=1.6605e-27
I0=current
pi=3.14159265359
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   
freq0=beta0*clight/C

; Oeffnen der Datenfiles:

num_p=long(0)
openr, 1, data_dir+"patric.dat"
tem=fltarr(18)   
while(not EOF(1)) do begin
 readf,1,tem
 num_p=num_p+1
end
close,1

print, "Ncells:",num_p/Nelements

if (nmax*Nelements LT num_p AND nmax GT 0) then num_p=Nelements*nmax

openr, 1, data_dir+"patric.dat"
ppic=fltarr(18,num_p) 
readf,1,ppic
close,1

; Grid:

s=fltarr(num_p)
s(0:num_p-1)=ppic(0,0:num_p-1)

; offsets

offsetx=fltarr(num_p)
offsety=fltarr(num_p)
pickupx=fltarr(num_p)
pickupy=fltarr(num_p)
xz=fltarr(num_p)

offsetx(0:num_p-1)=ppic(13,0:num_p-1)
offsety(0:num_p-1)=ppic(14,0:num_p-1)
pickupx(0:num_p-1)=ppic(16,0:num_p-1) 
pickupy(0:num_p-1)=ppic(17,0:num_p-1)
xz(0:num_p-1)=ppic(9,0:num_p-1)

num_win=floor(num_p/nwin)
offsetx_win=fltarr(nwin,num_win)
pickupx_win=fltarr(nwin,num_win)
for j=0, nwin-1 do begin
   for l=0L, num_win-1 do begin
       offsetx_win[j,l]=offsetx[j*num_win+l]
       pickupx_win[j,l]=pickupx[j*num_win+l]
   endfor
endfor   

; FFT of offsets:

han=hanning(num_p,/double)

offsetx_fft=fft(han*offsetx,/double)
offsety_fft=fft(han*offsety,/double)
pickupx_fft=fft(han*pickupx,/double)
pickupy_fft=fft(han*pickupy,/double)
xz_fft=fft(han*xz,/double)

Dt=(s(1)-s(0))/(beta0*clight)  ;  time step in track
freq=findgen(num_p/2+1) / (num_p*Dt)  ; resulting frequency grid
tunevec=freq/(beta0*clight/C)

han_win=hanning(num_win,/double)

offsetx_win_fft=dcomplexarr(nwin,num_win)
pickupx_win_fft=dcomplexarr(nwin,num_win)

for j=0, nwin-1 do begin
    offsetx_win_fft[j,*]=fft(han_win*offsetx_win[j,*],/double)
    pickupx_win_fft[j,*]=fft(han_win*offsetx_win[j,*],/double)
endfor    

freq_win=findgen(num_win/2+1) / (num_win*Dt)  ; resulting frequency grid
tunevec_win=freq_win/(beta0*clight/C)

offsetx_avg_fft=offsetx_win_fft[0,*]/nwin
pickupx_avg_fft=pickupx_win_fft[0,*]/nwin
for j=1, nwin-1 do begin
  offsetx_avg_fft[*]=offsetx_avg_fft[*]+offsetx_win_fft[j,*]/nwin
  pickupx_avg_fft[*]=pickupx_avg_fft[*]+pickupx_win_fft[j,*]/nwin
endfor

; dipole excitation signal:

dtheta=fltarr(num_p)
dtheta(0:num_p-1)=ppic(15,0:num_p-1)
dtheta_fft=fft(han*dtheta,/double)

; plot offset:

if( offset eq 1) then begin

maxoffset=max(offsetx(0:num_p-1))
baretune=Ncell*97.3/360.0

plot, s/(beta0*clight)*1.0e3, offsetx(0:num_p-1),/xstyle,/ystyle,xtitle="t [ms]",ytitle="offset in x"
;oplot, s/(beta0*clight)*1.0e3, offsety(0:num_p-1), color=150
oplot,s/(beta0*clight)*1.0e3, pickupx(0:num_p-1), line=2

endif

if( fft_offset eq 1) then begin

nmode=0.0
baretune=Ncell*97.3/360.0
maxoffset=max(smooth(abs(offsetx_fft(0:num_p/2))^2,1))
maxxz=max(smooth(abs(xz_fft(0:num_p/2))^2,1))
maxtheta=max(abs(dtheta_fft(0:num_p/2))^2)
;maxoffset=max(abs(offsetx_fft(0:num_p/2))^2)
maxtune=freq[!C]/(beta0*clight/C)
minoffset=5.0e-3*maxoffset
Seff=abs((0.5*baretune+nmode)*sqrt(5.0)*6.0e-3)  ;    0.009 

minrel=1.0e-3

plot,tunevec(0:num_p/2),smooth(abs(pickupx_fft(0:num_p/2))^2,1),/ylog,$
/xstyle,/ystyle,title="",xtitle="f/f!L0!N",ytitle="P [arb. units]",xrange=[nmode+baretune-2.0*abs(Seff),nmode+baretune+2.0*abs(Seff)],$xrange=[1.59,1.65],$xrange=[3.2,3.3],$
yrange=[minrel*maxoffset,maxoffset*2.0],xmargin=[12,3];,/nodata

;oplot, tunevec(0:num_p/2), maxoffset/maxtheta*abs(dtheta_fft(0:num_p/2))^2,color=100

;legend,['PATRIC','P(f)','P!U0!N(f)','Q!L0!N+!7D!3Q!Usc!N'],linestyle=[0,1,2,1],colors=[0,100,200,50],$
;box=0,charsize=1.2,/right,spacing=1.5

;oplot,tunevec(0:num_p/2),smooth(abs(pickupy_fft(0:num_p/2))^2,3),color=100,line=2

;oplot,tunevec(0:num_p/2),smooth(abs(offsetx_fft(0:num_p/2))^2,1),color=100
;oplot,tunevec(0:num_p/2),maxoffset/maxxz*smooth(abs(xz_fft(0:num_p/2))^2,1),color=200

;oplot,tunevec_win(0:num_win/2),smooth(abs(offsetx_avg_fft(0:num_win/2))^2,3),color=150

oplot,tunevec(0:num_p/2),maxoffset*(exp(-0.5*(tunevec-baretune)^2/(Seff*0.5)^2)),line=2,color=100

cohtune=10.0-2.0*baretune
incohtune=nmode+baretune  
syntune=8.2e3/freq0 ;   8.2e3/freq0 ;   8.0

print, "syntune:", syntune
print, baretune
print, maxtune

dQsc=0.0

fdist_vec=tunevec(0:num_p/2)
for j=0L, num_p/2 do fdist_vec[j]=fdist(tunevec[j]-incohtune,Seff)

plotS, [incohtune,incohtune],[minrel*maxoffset,2.0*maxoffset],line=1;,color=100
;plotS, [incohtune-syntune-dQsc,incohtune-syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+syntune-dQsc,incohtune+syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune-2*syntune-dQsc,incohtune-2*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+2*syntune-dQsc,incohtune+2*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune-3*syntune-dQsc,incohtune-3*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+3*syntune-dQsc,incohtune+3*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune-4*syntune-dQsc,incohtune-4*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+4*syntune-dQsc,incohtune+4*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune-5*syntune-dQsc,incohtune-5*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+5*syntune-dQsc,incohtune+5*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune-6*syntune-dQsc,incohtune-6*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+6*syntune-dQsc,incohtune+6*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune-7*syntune-dQsc,incohtune-7*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200
;plotS, [incohtune+7*syntune-dQsc,incohtune+7*syntune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2,color=200


for j=1, 10 do begin
  plotS, [j-incohtune,j-incohtune],[minrel*maxoffset,2.0*maxoffset],line=1,color=200
  plotS, [j+incohtune,j+incohtune],[minrel*maxoffset,2.0*maxoffset],line=1
endfor
  
plotS, [-1.0-incohtune,-1.0-incohtune],[minrel*maxoffset,2.0*maxoffset],line=1,color=200
plotS, [-1.0+incohtune,-1.0+incohtune],[minrel*maxoffset,2.0*maxoffset],line=1

dQsc=0.009/2.0*0.0

;oplot,tunevec(0:num_p/2)-dQsc,maxoffset*fdist_vec, line=3, color=200

; Blas. model



plotS, [incohtune-dQsc,incohtune-dQsc],[minrel*maxoffset,2.0*maxoffset],line=2

plotS, [incohtune+Qshift(dQsc,syntune,1),incohtune+Qshift(dQsc,syntune,1)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,-1),incohtune+Qshift(dQsc,syntune,-1)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,2),incohtune+Qshift(dQsc,syntune,2)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,-2),incohtune+Qshift(dQsc,syntune,-2)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,3),incohtune+Qshift(dQsc,syntune,3)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,-3),incohtune+Qshift(dQsc,syntune,-3)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,4),incohtune+Qshift(dQsc,syntune,4)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100
plotS, [incohtune+Qshift(dQsc,syntune,-4),incohtune+Qshift(dQsc,syntune,-4)],[minrel*maxoffset,2.0*maxoffset],line=2,color=100

endif

if(druck EQ 1) then begin
 device,/close
 spawn,'display schottky.eps' 
 set_plot,"x"
endif

end
