read_sdo,'hmi.sharp_cea_720s.2040.20120928_001200_TAI.Br.fits',index,data
s=size(data)
mpix_size=0.5;
nx=s(1)
ny=s(2)
bz=data
lylxbz=shift(double(bz),-1,-1)
lxbz=shift(double(bz),-1,0)
lybz=shift(double(bz),0,-1)
dBzdx = 0.5d0*(lylxbz+lxbz-lybz-double(bz))
dBzdy = 0.5d0*(lylxbz+lybz-lxbz-double(bz))

 arcsec2m = 696.0e3 
 units = 1.0/(arcsec2m*mpix_size)  ; ~ 0.0193*5.656/pixsize

   fdBzdx=dBzdx*units
   fdBzdy= dBzdy*units
; Set the edges to zero since the above shifts will give garbage at edges

      fdBzdX(nx-1,*)=0 & fdBzdX(*,ny-1)=0
      fdBzdY(nx-1,*)=0 & fdBzdY(*,ny-1)=0

    sbxx= fdBzdX*fdBzdX
     sbyy= fdBzdY*fdBzdY
   gbz=sqrt(sbxx+sbyy)

    n=0L
    cutof= 50e-6  ; units are in Gauss/meter
   for i=0L, (nx-1) do begin
   for j= 0L, (ny-1) do begin

   if gbz(i,j) gt cutof then begin
      n=n+1
   endif

   end
   end

print,n
; Length of strong gradient neutral line.

  LGNL= n*mpix_size*arcsec2m
  
  print, LGNL
 
end


