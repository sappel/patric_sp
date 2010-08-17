;+-----------------------------------------------------------------------
; NAME:
;       findpeaks
; PURPOSE:
;       Find peaks in a data array. To be a peak a local maximum must have
;	a 25% higher range than the first and last quarter of the peak range.
; CATEGORY:
;       Data Analysis
; CALLING SEQUENCE:
;	IDL> peaks = findpeaks( data )
; INPUTS:
;	data - array 
; OUTPUTS:
;	peaks - (up to "maxPeaks") values of peaks found. Returns 0 if no
;		peaks found.
; KEYWORDS:
;   Optional Input:
;	maxPeaks - max # of peaks to search for
;	AboveFrac - Fraction of peak height that peak must be above the 
;		peak "tails" to be counted (default=0.25)
;	NeighborhoodFrac - fraction of whole X range to be a neighborhood.
;          	Range is initially divided up in these neighborhoods
;		for determination of initial maxima.
;	Plot - if set, will plot data with vertical lines indicating which
;		peaks qualified.
;	Ordered - if set, returns peaks in order of highest first
;   Optionally Returned:
;	indicesOfPeaks - array indices of peaks returned (=-1 if no peaks found)
;	FWHM - full width at half max estimate, in indices
; EXAMPLE:
;	IDL> data = mk_2peaks()
;	IDL> peaks = findpeaks( data, /plot )
; COMMON BLOCKS:
;       none
; NOTES:
;	Maxima too close to the beginning and end of array are not counted.
;       You may wish to do filtering or smoothing before passing data.
; LIMITATIONS:
;	What constitutes a "peak" is subjective, so this won't fit all needs
; MODIFICATION HISTORY:
;	08-May-2003 Written by Bill Davis, PPPL
;------------------------------------------------------------------------
function findpeaks, data, maxPeaks=maxPeaks, Plot=Plot,   $
                    AboveFrac=AboveFrac, NeighborhoodFrac=NeighborhoodFrac,  $
		    indicesOfPeaks=indicesOfPeaks,  $
		    FWHM=FWHM, ORDERED=Ordered, status=status

if n_elements(maxPeaks) eq 0 then maxPeaks=2
if n_elements(AboveFrac) eq 0 then AboveFrac = 0.25

npts = n_elements( data )

if n_elements(NeighborhoodFrac) eq 0 then begin
   if npts lt 100 then NeighborhoodFrac = 0.2 $
                  else NeighborhoodFrac = 0.1
endif

status = 0

nPtsRange = (npts * NeighborhoodFrac)>20
rangeInds = lindgen( nPtsRange )
bandInds = lindgen( nPtsRange/4 )
nRanges = NINT( npts/nPtsRange )

iStarts = NINT( LINDGEN( nRanges ) * nPtsRange)

maxLocs = LONARR( nRanges )
maxVals = FLTARR( nRanges )
minVals = FLTARR( nRanges )

for i=0, nRanges-1 do begin
   maxVals[i] = MAX( data[ rangeInds + iStarts[i] ],  iLoc, MIN=minVal )
   minVals[i] = minVal
   maxLocs[i] = iLoc+ iStarts[i]
endfor

  ; remove peaks too close to edge

leftOffset = nint(3./4.*nPtsRange/2)
rightOffset = nint(1./4.*nPtsRange/2)
OKInds = where( maxLocs-rightOffSet ge 0 and  $
                 maxLocs+rightOffSet le npts-1, count)
if count le 0 then begin
   indicesOfPeaks = -1
   return, [0.0]
endif

maxLocs = maxLocs[ OKInds ]
goodPeak = intarr( n_elements(maxLocs) )
 
bandInds = LINDGEN( NINT(nPtsRange/3) )

   ; first, make sure it's the highest in the immediate area
      
for i=0, n_elements(maxLocs)-1 do begin

   ImmedRange =  bandInds-(n_elements(bandInds)/2) +  maxLocs[i]
   MaxImmed = MAX( data[ 0>ImmedRange<(npts-1) ] )
   if maxImmed eq data[ maxLocs[i] ] then goodPeak[i]=1

endfor

   ; see if really a peak, i.e., is significantly above a base

goodInds = where( goodPeak eq 1, count )
if count le 0 then begin
   indicesOfPeaks = -1
   return, [0.0]
endif

maxLocs = maxLocs(goodInds)

goodPeak = lonarr( n_elements(maxLocs) )
FWHM = lonarr( n_elements(maxLocs) )

for i=0, n_elements(maxLocs)-1 do begin

   if i eq 0 then begin
      leftOffset = (2.0/3.0*maxLocs[i] )
      if n_elements(maxLocs) eq 1   $
         then rightOffset = maxLocs[i] + 2./3*(npts-maxLocs[i])  $
         else rightOffset = (maxLocs[i+1]-maxLocs[i])/2
   endif else if i eq n_elements(maxLocs)-1 then begin
      leftOffset = (maxLocs[i]-maxLocs[i-1])/2
      rightOffset = 2.0/3.0*( npts - maxLocs[i] )
   endif else begin
      leftOffset = (maxLocs[i]-maxLocs[i-1])/2
      rightOffset = (maxLocs[i+1]-maxLocs[i])/2
   endelse
   leftOffset = NINT(leftOffset)
   rightOffset = NINT(rightOffset)
   
      ; find minimum for heavily smoothed right & left regions
   nSmRight = (rightOffset/5) > 3
   rtRange = indgen(rightOffset) + NINT(maxLocs[i])
   if rightOffset ge 3 then begin
      smData = smooth(data[rtRange],nSmRight )
   endif else smData = data[rtRange]
   RtMin = MIN( smData, iRtLoc)
   rightOffset = rightOffset < (iRtLoc+1)

   nSmLeft = (leftOffset/5) > 3
   leftRange = NINT(maxLocs[i] - (leftOffset-indgen(leftOffset)))
   if leftOffset ge 3 then begin
      smData = smooth(data[leftRange],nSmLeft )
   endif else smData = data[leftRange]
   LeftMin = MIN( smData, iLeftLoc)
   leftOffset = leftOffset - iLeftLoc

      ; just consider range between local minima
   thisRange = indgen(rightOffset+leftOffset+1)
   thisRange = thisRange + maxLocs[i] - leftOffset 
   
      ;(if minima were both at ends, no need for the following)
   if iLeftLoc eq 0 and iRtLoc eq n_elements(rtRange)-1 then begin
      goodPeak[i] = 1
   endif else begin
      MaxRange = MAX( data[thisRange], MIN=MinRange )

      leftBandInds = (bandInds +  maxLocs[i] - leftoffset) > 0
      maxLeftBand = MAX( data[leftBandInds], MIN=minLeftBand )
      LeftHeight = (maxLeftBand-MinRange)

      rightBandInds = (bandInds +  maxLocs[i] + rightOffset) < (npts-1)
      maxRightBand = MAX( data[RightBandInds], MIN=minRightBand )
      RightHeight = (maxRightBand-MinRange)

      PeakHeight = (MaxRange-MinRange)
      if PeakHeight gt (1+AboveFrac)*LeftHeight AND $
	 PeakHeight gt (1+AboveFrac)*RightHeight  THEN goodPeak[i] = 1
   endelse
   
   FWHM[i] = max( thisRange ) - min( thisRange )
   
endfor

;;;stop

goodInds = where( goodPeak eq 1, count )
if count le 0 then begin
   indicesOfPeaks = -1
   return, [0.0]
endif

maxLocs = maxLocs[ goodInds ]
FWHM = FWHM[ goodInds ]

if keyword_set( plot ) then begin
    
   plot, data
   for iline=0,n_elements(maxlocs)-1 do $
       oplot,[maxlocs[iline],maxlocs[iline]],!y.crange, color=mk_color('green')
endif

   ; find highest peaks, up to the max number requested

IF n_elements(maxlocs) eq 0 THEN BEGIN
   indicesOfPeaks = -1
   return, [0.0]
endif

maxes = data[ maxLocs ]
isort = REVERSE( sort(maxes) )
sortedMaxes = maxes( iSort )

if maxPeaks le n_elements( maxLocs ) then begin
   peaks = sortedMaxes[ 0:maxPeaks-1 ]
   indicesOfPeaks = maxLocs[ iSort[0:maxPeaks-1] ]
   FWHM = FWHM[ iSort[0:maxPeaks-1] ]
endif else begin
   peaks = sortedMaxes
   indicesOfPeaks = maxLocs  
endelse

if NOT KEYWORD_SET( Ordered ) then begin
   iSort = sort( indicesOfPeaks )
   indicesOfPeaks = indicesOfPeaks[iSort]
   peaks = peaks[iSort]
endif

status=1
return, peaks
end
