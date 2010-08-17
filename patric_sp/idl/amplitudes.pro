pro amplitudes, data_dir

!P.MULTI=[0,1,1]
!P.CHARSIZE=1.7
!P.CHARTHICK=3.0
!P.THICK=5.0
!Y.THICK=5.0
!X.THICK=5.0

;Eingabeparameter:

pict_max=301                    ; Max. Anzahl an Bildern
movie=0                         ; movie (gif/screen/nein) (2/1/0)
druck=1                         ; Ausdruck ja/nein (1/0)
loadct,39                       ; Color: 39 (Rainbow), 3 (red temperature)
dpmax=1.0e-2                    ; dp/p max
xmax=0.01                       ; x max (m)
ymax=0.01                       ; y max (m)
ds=0.4                          ; determines xsmax and ysmax
smoothing=6                     ; contour smoothing
clevels=40                      ; contour levels
hmax=10

set_plot,"x"
if(druck EQ 1) then begin
  set_plot,"ps"
  device,filename="amplitudes.eps",/color,/encapsulated,/inches,ysize=4.8 ; 9.6
endif 
;device,decomposed=0
;x_size=640
;y_size=512

; read parameter file idl.dat:

openr, 1, data_dir+"idl.dat"
readf,1,numprocs                ; Prozesse
readf,1,e_kin                   ; kinetische Energie/MeV
readf,1,qb                      ; Ladung/qe
readf,1,mb                      ; Masse/mp
readf,1,current                 ; Strom      
readf,1,C                       ; Ringumfang
readf,1,Nelements               ; elements in cell 
readf,1,cell_length             ; Zellenlaenge
readf,1,Np                      ; Teilchen fuer Plot
readf,1,pipe_radius
readf,1,NX
readf,1,NY
readf,1,NZ
close,1

; Constants:

qe=1.6022e-19
mp=1.6605e-27
I0=current
pi=3.14159265359
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))
Dt = cell_length/(beta0*clight)

; get num_p

pics=fltarr(pict_max,8,Np*numprocs) 

num_p=0L
openr, 1, data_dir+"pics_0.dat"
tem=fltarr(8,Np)   
while(not EOF(1)) do begin
  readu,1,tem
  pics(num_p,0:7,0:Np-1)=tem(0:7,0:Np-1)
  num_p=num_p+1
end
close,1

; dipole array 

dipole_x=fltarr(num_p,NZ)

; get dipole data

openr, 3, data_dir+"dipole_x.dat"
tem=fltarr(NZ)
for i=0,num_p-1 do begin
  readu,3,tem
  dipole_x(i,0:NZ-1)=tem
endfor
close,3

; FFT

print_cell=120.0

dipole_x[0,*]=dipole_x[1,*]

dipolex_fft=fltarr(num_p,NZ/2+1)
h=findgen(NZ/2+1)
t=findgen(num_p)
t=t*print_cell*Dt

for j=0, num_p-1 do begin
  u=FFT(dipole_x[j,0:NZ-1]/I0)
  dipolex_fft[j,0:NZ/2]=2.0*ABS(u[0:NZ/2])
endfor

; plot

;plot, t*1.0e3, dipolex_fft[*,3],/xstyle,/ystyle,xtitle="t [ms]"
;oplot, t*1.0e3, 0.5*dipolex_fft[0,3]*exp(t*1.0e3/0.023), line=1

clevels=40
dmax=max(dipolex_fft[*,*])
dmin=min(dipolex_fft[*,*])
contour,dipolex_fft[*,*],t*1.0e3,h,/xstyle,/ystyle,$
        xtitle="t [ms]",ytitle="harmonic",$
        nlevels=clevels,/fill,$
        levels=dmin+findgen(clevels+1)*(dmax-dmin)/clevels,$
        C_COLORS=25+findgen(clevels+1)*(255-25)/clevels,yrange=[0,hmax]

if(druck EQ 1) then begin
  device,/close
  spawn,'display amplitudes.eps' 
  set_plot,"x"
endif

END

