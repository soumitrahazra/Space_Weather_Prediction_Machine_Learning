;+
;NAME:
;     get_neutral
;PURPOSE:
;     Computes a smoothed field suitable for contouring a neutral line from an 
;     input longitudinal magnetic field array.
;CATEGORY:
;CALLING SEQUENCE:
;     neutral = get_neutral(B_long)
;INPUTS:
;     Array of longitudinal magnetic field
;OPTIONAL INPUT PARAMETERS:
;     smooth = degree of smoothing, the bigger the number, the smoother the
;              neutral line  (default=0.25; 0.0=no smoothing)
;KEYWORD PARAMETERS
;OUTPUTS:
;     Array of values suitable for contouring a neutral line
;COMMON BLOCKS:
;SIDE EFFECTS:
;RESTRICTIONS:
;     The neutral line is smoothed and hence IS NOT EXACT.
;PROCEDURE:
;  Example
;   IDL>  n = get_neutral(b.b_long)
;   IDL>  nx = n_elements(n(*,0))
;   IDL>  contour,n,levels=[0.0],c_line=[0],c_labels=[0],spline=1.0/nx,c_thick=5
;MODIFICATION HISTORY:
;     TRM 4/1992
;     TRM 12/1992 Added smooth parameter to control the amount of smoothing
;-

function get_neutral,b,smooth,grow=grow

  if NOT keyword_set(grow) then begin
     if n_elements(smooth) LE 0 then smooth = 0.25d0
     nx = n_elements(b(*,0))
     ny = n_elements(b(0,*))
     neutral = double(b)
     neutral(0,*) = 0      ; Since shift wraps around, set the edges to zero
     neutral(nx-1,*) = 0
     neutral(*,0) = 0
     neutral(*,ny-1) = 0
     for i=-1,1 do begin
        for j=-1,1 do begin
           if i NE 0 OR j NE 0 then neutral=neutral + shift(b,i,j)*smooth
        endfor
     endfor

     return, neutral
  endif $
  else begin

     ; Binary grow algorithm

     nx = n_elements(b(*,0))
     ny = n_elements(b(0,*))
     neutral = fltarr(nx+2,ny+2)
     neutral(1:nx,1:ny) = float(b)
  
     bmax = shift(neutral,1,0) > $
           shift(neutral,0,1)  > $
           shift(neutral,-1,0) > $
           shift(neutral,0,-1) > $
           shift(neutral,1,1) > $
           shift(neutral,-1,1) > $
           shift(neutral,1,-1) > $
           shift(neutral,-1,-1) > $
           neutral

     bmin = shift(neutral,1,0) < $
           shift(neutral,0,1)  < $
           shift(neutral,-1,0) < $
           shift(neutral,0,-1) < $
           shift(neutral,1,1) < $
           shift(neutral,-1,1) < $
           shift(neutral,1,-1) < $
           shift(neutral,-1,-1) < $
           neutral

     neutral = neutral + bmin + bmax
     return,neutral(1:nx,1:ny)
  endelse

end
