pro tphasespace, data_dir

!P.MULTI=[0,2,2]
!P.CHARSIZE=1.4
!P.CHARTHICK=3.0
!P.THICK=5.0
!Y.THICK=5.0
!X.THICK=5.0

;Eingabeparameter:

pict_max=301                    ; Max. Anzahl an Bildern
movie=0                         ; movie (gif/screen/nein) (2/1/0)
druck=0                         ; Ausdruck ja/nein (1/0)
loadct,39                       ; Color: 39 (Rainbow), 3 (red temperature)
dpmax=1.0e-2                    ; dp/p max
xmax=0.045                      ; x max (m)
ymax=0.045                      ; y max (m)
ds=0.4                          ; determines xsmax and ysmax
smoothing=3                     ; contour smoothing
clevels=40                      ; contour levels

set_plot,"x"
device,decomposed=0
x_size=640
y_size=512

if(druck EQ 1) then begin
  set_plot,"ps"
  device,filename="ppic.eps",/encapsulate,/inches,ysize=6.0,xsize=7.0,/color ;,bits=16
endif 

if(movie eq 0) then set_plot,"ps"
if(movie eq 0) then spawn,'rm -rf /tmp/lobo*'
if(movie eq 2) then spawn,'rm -rf /tmp/lobo*'
if(movie eq 2) then device,decomposed=1

; plot options (1/0), four must be chosen:

pz=0
rhoxy=1
rhoxxs=0	
rhoyys=0
rhoxsys=0
rhozx=1
phaseadvance=0
efieldxy=0
rhox=0
rhoz=0
dipolex=1

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

; fields and vectors:

pics=fltarr(pict_max,8,Np*numprocs) 

; Constants:

qe=1.6022e-19
mp=1.6605e-27
I0=current
pi=3.14159265359
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   

; Grids:

z=fltarr(NZ)
z(0)=-0.5*C
for i=1, NZ-1 do z(i)=z(i-1)+C/(NZ-1.0)

x=fltarr(NX)
x(0)=-pipe_radius
for i=1, NX-1 do x(i)=x(i-1)+2.0*pipe_radius/(NX-1.0)

xsmax=1.0e3*4.0*pipe_radius/(ds*(NX-1.0))
xs=fltarr(NX)
xs(0)=-xsmax
for i=1, NX-1 do xs(i)=xs(i-1)+2.0*xsmax/(NX-1.0)

y=fltarr(NY)
y(0)=-pipe_radius
for i=1, NY-1 do y(i)=y(i-1)+2.0*pipe_radius/(NY-1.0)

ysmax=1.0e3*4.0*pipe_radius/(ds*(NY-1.0))
ys=fltarr(NX)
ys(0)=-ysmax
for i=1, NY-1 do ys(i)=ys(i-1)+2.0*ysmax/(NY-1.0)


; Oeffnen der PIC Datenfiles:

num_p=0
openr, 1, data_dir+"pics_0.dat"
tem=fltarr(8,Np)   
while(not EOF(1)) do begin
  readu,1,tem
  pics(num_p,0:7,0:Np-1)=tem(0:7,0:Np-1)
  num_p=num_p+1
end
close,1

for i=1,numprocs-1 do begin
  openr, 1, data_dir+"pics_"+string(strtrim(i,1))+".dat"   
  for j=0, num_p-1 do begin
    readu,1,tem
    pics(j,0:7,i*Np:(i+1)*Np-1)=tem(0:7,0:Np-1)
  endfor
  close,1 
endfor

; define fields for momenta

rho_z=fltarr(num_p,NZ)
dipole_x=fltarr(num_p,NZ)
rho_xy=fltarr(num_p,NX,NY)
rho_xxs=fltarr(num_p,NX,NX)
rho_yys=fltarr(num_p,NY,NY)
rho_zx=fltarr(num_p,NZ,NX)
rho_xsys=fltarr(num_p,NX,NY)
Ex=fltarr(num_p,NX,NY)
Ey=fltarr(num_p,NX,NY)
rho_x=fltarr(num_p,NX)
rho_y=fltarr(num_p,NY)

; opening of momenta files

;if( rhoz eq 1 ) then begin
openr, 3, data_dir+"rho_z.dat"
tem=fltarr(NZ)
for i=0,num_p-1 do begin
  readu,3,tem
  rho_z(i,0:NZ-1)=tem
endfor
close,3
;endif


if( dipolex eq 1 ) then begin
  openr, 3, data_dir+"dipole_x.dat"
  tem=fltarr(NZ)
  for i=0,num_p-1 do begin
    readu,3,tem
    dipole_x(i,0:NZ-1)=tem
  endfor
  close,3
endif


if( rhoxxs eq 1 ) then begin
  openr, 3, data_dir+"xxs.dat"
  tem=fltarr(NX)
  for i=0,num_p-1 do begin
    for k=0,NX-1 do begin
      readu,3,tem
      rho_xxs(i,k,0:NX-1)=tem
    endfor
  endfor
  close,3
endif

if( rhoyys eq 1 ) then begin
  openr, 3, data_dir+"yys.dat"
  tem=fltarr(NY)
  for i=0,num_p-1 do begin
    for k=0,NY-1 do begin
      readu,3,tem
      rho_yys(i,k,0:NY-1)=tem
    endfor
  endfor
  close,3
endif

if( rhoxsys eq 1 ) then begin
  openr, 3, data_dir+"xsys.dat"
  tem=fltarr(NY)
  for i=0,num_p-1 do begin
    for k=0,NX-1 do begin
      readu,3,tem
      rho_xsys(i,k,0:NY-1)=tem
    endfor
  endfor
  close,3
endif

if( rhoxy eq 1 or rhox eq 1 ) then begin
  openr, 3, data_dir+"rho_xy.dat"
  tem=fltarr(NY)
  for i=0,num_p-1 do begin
    for k=0,NX-1 do begin
      readu,3,tem
      rho_xy(i,k,0:NY-1)=tem
    endfor
  endfor
  close,3
endif

if( rhozx eq 1 ) then begin
  openr, 3, data_dir+"zx.dat"
  tem=fltarr(NX)
  for i=0,num_p-1 do begin
    for k=0,NZ-1 do begin
      readu,3,tem
      rho_zx(i,k,0:NX-1)=tem
    endfor
  endfor
  close,3
endif

if( efieldxy eq 1 ) then begin
  openr, 3, data_dir+"Ex.dat"
  tem=fltarr(NY)
  for i=0,num_p-1 do begin
    for k=0,NX-1 do begin
      readu,3,tem
      Ex(i,k,0:NY-1)=tem
    endfor
  endfor
  close,3

  openr, 3, data_dir+"Ey.dat"
  tem=fltarr(NY)
  for i=0,num_p-1 do begin
    for k=0,NX-1 do begin
      readu,3,tem
      Ey(i,k,0:NY-1)=tem
    endfor
  endfor
  close,3
endif

pipe_radius=xmax                ; !!!!!!!!  overwrite pipe radius
xsmax=1.0*xsmax
ysmax=1.0*ysmax


;IDL device commands:

if (movie EQ 1) then begin
  xinteranimate,set=[x_size,y_size,num_p] ; set up animation
endif

;window,0,xsize=x_size,ysize=y_size

; start plot loop:

for j=0,num_p-1 do begin

; print eps's:
  
  if( movie eq 0) then begin
    if (j lt 9) then begin
      device,filename=strcompress("/tmp/lobo00" + string(j+1) + ".eps", $
                                  /remove_all),/color,/encapsulated
    endif
    if (j ge 9 AND j lt 99) then begin
      device,filename=strcompress("/tmp/lobo0" + string(j+1) + ".eps", $
                                  /remove_all),/color,/encapsulated
    endif
    if (j ge 99) then begin
      device,filename=strcompress("/tmp/lobo" + string(j+1) + ".eps", $
                                  /remove_all),/color,/encapsulated
    endif
  endif


;!p.title=strcompress("turn="+string(j),/remove_all)
  !p.title=strcompress("time [ms] = "+string(1.0e3*j*120.0*cell_length/(beta0*clight)),/remove_all)

  if(rhoz EQ 1) then begin
    plot,z,beta0*clight*rho_z[j,*]/I0,/xstyle,/ystyle,yrange=[0,10]
  endif

  if(dipolex EQ 1) then begin
    plot,z,dipole_x[j,*]/I0,/xstyle,/ystyle,yrange=[-0.25*xmax,0.25*xmax],ymargin=[3.5,2],$
         xmargin=[10,-8],ytitle="dipole current [a.u.]",XTICKFORMAT="(A1)" 
;oplot,z,0.2*xmax*beta0*clight*rho_z[j,*]/I0,line=1
    plot,z,dipole_x[j,*]/I0,yrange=[-0.5*xmax,0.5*xmax],title="",/nodata,xstyle=4,ystyle=4, $
         color=0
  endif

  if(pz EQ 1) then begin
    plot,pics(j,5,0:numprocs*Np-1),pics(j,4,0:numprocs*Np-1),$
         /xstyle,/ystyle,psym=3,xtitle="z",ytitle="dp/p",$
         xrange=[-0.5*C,0.5*C],yrange=[-dpmax,dpmax],/nodata 

    for i=0,numprocs-1 do begin
      oplot,pics(j,5,i*Np:(i+1)*Np-1),pics(j,4,i*Np:(i+1)*Np-1),psym=3,color=(i+1)*250/numprocs
    endfor
  endif

  if(rhozx EQ 1) then begin
    plot,pics(j,5,0:numprocs*Np-1),100.0*pics(j,0,0:numprocs*Np-1),$
         /xstyle,/ystyle,psym=2,symsize=0.05,xtitle="z [m]",ytitle="x [cm]",title=" ",$
         xrange=[-0.5*C,0.5*C],yrange=[-xmax*100,xmax*100],xmargin=[10,-8],ymargin=[4,-3.5], $
         /nodata

    for i=0,numprocs-1 do begin
      oplot,pics(j,5,i*Np:(i+1)*Np-1),100.0*pics(j,0,i*Np:(i+1)*Np-1),psym=2,symsize=0.05 ;,color=(i+1)*250/numprocs
    endfor
    arr=fltarr(NZ,NX)
    arr(*,*)=rho_zx(j,*,*)
    rho_zx(j,*,*)=smooth(arr,smoothing,/edge_truncate)
    picmax=max(abs(rho_zx(j,*,*)))
    picmin=picmax*1.0e-1
    contour,rho_zx(j,*,*),z,x*100.0,/xstyle,/ystyle,$
            xtitle="z [cm]",ytitle="x [cm]",title=" ",yrange=[-xmax*100,xmax*100],$
            xrange=[-0.5*C,0.5*C],nlevels=clevels,/fill,$
            levels=picmin+findgen(clevels+1)*(picmax-picmin)/clevels,xmargin=[10,-8], $
            ymargin=[4,-3.5], C_COLORS=25+findgen(clevels+1)*(255-25)/clevels,/overplot
  endif

  if(rhoxy EQ 1) then begin
    plot,pics(j,2,0:numprocs*Np-1)*100,pics(j,0,0:numprocs*Np-1)*100,$
         /xstyle,psym=2,symsize=0.05,$
         yrange=[-ymax*100,ymax*100],$
         xrange=[-pipe_radius*100,pipe_radius*100],title=" ",$
         xtitle="y [cm]",ytitle="x [cm]",xmargin=[8,4.5],ystyle=5,ymargin=[4,-3.5] ;,xticks=3
    axis,yaxis=1,/ystyle
    arr=fltarr(NX,NY)           ; !!!!!!   NX=NY !!!!!!
    for i=0, NX-1 do begin
      for k=0, NX-1 do begin
        arr[i,k]=rho_xy[j,k,i]
      endfor
    endfor  
    rho_xy(j,*,*)=smooth(arr,smoothing,/edge_truncate)
    picmax=max(abs(rho_xy(j,*,*)))
    picmin=picmax*1.0e-1
    contour,rho_xy(j,*,*),y*100,x*100,/xstyle,/ystyle,title=" ",$
            xtitle="z [cm]",ytitle="x [cm]",yrange=[-ymax*100,ymax*100],$
            xrange=[-pipe_radius*100,pipe_radius*100],nlevels=clevels,/fill,$
            levels=picmin+findgen(clevels+1)*(picmax-picmin)/clevels,$
            C_COLORS=25+findgen(clevels+1)*(255-25)/clevels,/overplot
  endif

  if(rhox EQ 1) then begin
    rhox_vec=fltarr(NX)
    for l=0, NX-1 do begin
      rhox_vec[l]=0 
      for i=0, NY-1 do begin
        rhox_vec[l]+=rho_xy[j,l,i]
      endfor  
    endfor

    rhox_norm=0.0
    for l=0, NX-1 do begin
      rhox_norm+=rhox_vec[l]*2.0*pipe_radius/(NX-1.0)
    endfor

    plot,rhox_vec(*)/rhox_norm,x*100.0,/xlog,/xstyle,$
         xtitle="density [arb. units]",xrange=[0.02,90.0], $
         yrange=[-pipe_radius*100,pipe_radius*100],$
         xmargin=[8,4.5],ystyle=5,ymargin=[4,-3.5],title="",XtickFormat='(F10.1)'
    axis,yaxis=1,/ystyle
  endif

  if(rhoxxs EQ 1) then begin
    plot,pics(j,1,0:numprocs*Np-1)*1000.0,pics(j,0,0:numprocs*Np-1)*100.0,$
         /xstyle,psym=2,symsize=0.05,yrange=[-pipe_radius*100,pipe_radius*100],$
         xrange=[-xsmax,xsmax],xtitle="x' [mrad]",title="",xmargin=[8,4.5],ystyle=5,ymargin=[4,-3.5]
    axis,yaxis=1,/ystyle
    arr=fltarr(NX,NY)           ; !!!!!!   NX=NY !!!!!!
    for i=0, NX-1 do begin
      for k=0, NX-1 do begin
        arr[i,k]=rho_xxs[j,k,i]
      endfor
    endfor  
    rho_xxs(j,*,*)=smooth(arr,smoothing,/edge_truncate)
    picmax=max(abs(rho_xxs(j,*,*)))
    picmin=picmax*1.0e-1
    contour,rho_xxs(j,*,*),xs,x*100,/xstyle,/ystyle,$
            xtitle="x [cm]",ytitle="x' [mrad]",xrange=[-xsmax,xsmax],$
            yrange=[-pipe_radius*100,pipe_radius*100],nlevels=clevels,/fill,$
            levels=picmin+findgen(clevels+1)*(picmax-picmin)/clevels,$
            C_COLORS=20+findgen(clevels+1)*(254-20)/clevels,/overplot
  endif

  if(rhoyys EQ 1) then begin
    plot,pics(j,2,0:numprocs*Np-1)*100,pics(j,3,0:numprocs*Np-1)*1000.0,$
         /xstyle,/ystyle,psym=3,xrange=[-ymax*100,ymax*100],$
         yrange=[-ysmax,ysmax],xtitle="y [cm]",ytitle="y' [mrad]"
    arr=fltarr(NY,NY)
    arr(*,*)=rho_yys(j,*,*)
    rho_yys(j,*,*)=smooth(arr,smoothing,/edge_truncate)
    picmax=max(abs(rho_yys(j,*,*)))
    picmin=picmax*1.0e-1
    contour,rho_yys(j,*,*),y*100,ys,/xstyle,/ystyle,$
            xtitle="x [cm]",ytitle="x' [mrad]",yrange=[-ysmax,ysmax],$
            xrange=[-ymax*100,ymax*100],nlevels=clevels,/fill,$
            levels=picmin+findgen(clevels+1)*(picmax-picmin)/clevels,$
            C_COLORS=20+findgen(clevels+1)*(254-20)/clevels,/overplot
  endif

  if(rhoxsys EQ 1) then begin
    plot,100.0*pics(j,1,0:numprocs*Np-1),100.0*pics(j,3,0:numprocs*Np-1),$
         /xstyle,/ystyle,psym=3,xtitle="x' [mrad]",ytitle="y' [mrad]",$
         xrange=[-xsmax,xsmax],yrange=[-ysmax,ysmax] 
    arr=fltarr(NX,NY)
    arr(*,*)=rho_xsys(j,*,*)
    rho_xsys(j,*,*)=smooth(arr,smoothing,/edge_truncate)
    picmax=max(abs(rho_xsys(j,*,*)))
    picmin=picmax*1.0e-1
    contour,rho_xsys(j,*,*),xs,ys,/xstyle,/ystyle,$
            xtitle="z' [mrad]",ytitle="x' [mrad]",$
            yrange=[-ysmax,ysmax],xrange=[-xsmax,xsmax],nlevels=clevels,/fill,$
            levels=picmin+findgen(clevels+1)*(picmax-picmin)/clevels,$
            C_COLORS=25+findgen(clevels+1)*(254-25)/clevels,/overplot
  endif

  if( phaseadvance eq 1 ) then begin
    Qx0=5.0                     ;3.5
    Qx1=7.0                     ;4.5
    Qy0=5.0                     ;3.5
    Qy1=7.0                     ;4.5
    Qxbare=6.2                  ;4.23
    Qybare=6.2                  ;4.23

    Qycoh=Qybare
    Qxcoh=Qxbare-0.0

    C2=100.0*!PI*2.0

    plot,0.5*C2/pics(j,6,0:numprocs*Np-1),0.5*C2/pics(j,7,0:numprocs*Np-1),$
         /xstyle,/ystyle,psym=3,xtitle="tunex",ytitle="tuney",$
         xrange=[Qx0,Qx1],yrange=[Qy0,Qy1],/nodata
    for i=0, numprocs-1 do begin
      oplot,0.5*C2/pics(j,6,i*Np:(i+1)*Np-1),0.5*C2/pics(j,7,i*Np:(i+1)*Np-1), $
            psym=3,color=(i+1)*250/numprocs
    endfor
;plot,pics(j,6,0:numprocs*Np-1)*12.0/(2.0*!PI),pics(j,7,0:numprocs*Np-1)*12.0/(2.0*!PI),$
;/xstyle,/ystyle,psym=3,xtitle="tunex",ytitle="tuney",$
;xrange=[Qx0,Qx1],yrange=[Qy0,Qy1],/nodata
;for i=0, numprocs-1 do begin
;   oplot,pics(j,6,i*Np:(i+1)*Np-1)*12.0/(2.0*!PI),pics(j,7,i*Np:(i+1)*Np-1)*12.0/(2.0*!PI),psym=3,color=(i+1)*250/numprocs
;endfor
    plots,Qxbare,Qybare,psym=4
    plotS, [Qxcoh,Qxcoh],[Qy0,Qy1],line=2
  endif

  if( efieldxy eq 1) then begin
;plot,x*100,Ex(j,*,NY/2),/xstyle,xrange=[-pipe_radius*100,pipe_radius*100]
    contour,Ex[j,*,*]^2+Ey[j,*,*]^2,x*100,y*100,/xstyle,/ystyle,nlevels=clevels,/fill
  endif

; load animation:

  if (movie EQ 1) then begin
    xinteranimate,frame=j,window=0
  endif

; print gif's:

  if( movie eq 2) then begin
    if (j lt 9) then write_bmp,$
      strcompress("/tmp/lobo00" + string(j+1) + ".bmp",/remove_all),tvrd()
    if (j ge 9 AND j lt 99) then write_bmp,$
      strcompress("/tmp/lobo0" + string(j+1) + ".bmp",/remove_all),tvrd()
    if (j ge 99) then write_bmp,$ 
      strcompress("/tmp/lobo" + string(j+1) + ".bmp",/remove_all),tvrd()
  endif
endfor

if (movie EQ 1) then begin
  wdelete,0
endif


;play animation:

if (movie EQ 1) then begin
  xinteranimate
endif

if(druck EQ 1) then begin 
  device,/close
  spawn,'gv ppic.eps'
  set_plot,"x"
endif

end
