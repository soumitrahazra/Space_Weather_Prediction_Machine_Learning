"""
Purpose:   Calculation of solar magnetic parameters using SDO/HMI vector magnetic field data:

           USFLUX  Total unsigned flux in Maxwells
           CMASK   Number of pixels used in the USFLUX calculation
           MEANGAM Mean inclination angle, gamma, in degrees
           TOTBSQ  Total magnitude of Lorentz force
           TOTFX   Sum of X-component of Lorentz force
           TOTFY   Sum of Y-component of Lorentz force
           TOTFZ   Sum of Z-component of Lorentz force
           EPSX    Sum of X-component of normalized Lorentz force
           EPSY    Sum of Y-component of normalized Lorentz force
           EPSZ    Sum of Z-component of normalized Lorentz force 
           MEANGBT Mean value of the total field gradient, in Gauss/Mm
           MEANGBZ Mean value of the vertical field gradient, in Gauss/Mm
           MEANGBH Mean value of the horizontal field gradient, in Gauss/Mm
           MEANJZD Mean vertical current density, in mA/m2
           TOTUSJZ Total unsigned vertical current, in Amperes
           MEANALP Mean twist parameter, alpha, in 1/Mm
           MEANJZH Mean current helicity in G2/m
           TOTUSJH Total unsigned current helicity in G2/m
           ABSNJZH Absolute value of the net current helicity in G2/m
           SAVNCPP Sum of the absolute value of the net current per polarity in Amperes
           MEANPOT Mean photospheric excess magnetic energy density in ergs per cubic centimeter
           TOTPOT  Total photospheric magnetic energy density in ergs per centimeter
           MEANSHR Mean shear angle (measured using Btotal) in degrees
           SHRGT45 Area with shear angle greater than 45 degrees (as a percent of total area)
           R-VALUE Sum of flux near polarity inversion line
                      
Inputs:    We use the following segments:
           
           [example filename]                 --> [description]
           hmi.sharp_cea_*.Br.fits            --> radial component of the magnetic field vector
           hmi.sharp_cea_*.Bt.fits            --> theta-component of the magnetic field vector
           hmi.sharp_cea_*.Bp.fits            --> phi-component of the magnetic field vector
           hmi.sharp_cea_*.conf_disambig.fits --> bits indicate confidence levels in disambiguation result
           hmi.sharp_cea_*.bitmap.fits        --> bits indicate result of automatic detection algorithm
           
Usage:     This code depends on the numpy, scipy, and sunpy libraries.

Written:   Soumitra Hazra
Last Modified Date: August 26 2020
           
"""

# import some modules
import sunpy, sunpy.map, scipy, numpy as np, sys, math, argparse, astropy, glob
import pandas as pd
import os, os.path
from astropy.io import ascii
# define some constants
radsindeg = np.pi/180.
munaught  = 0.0000012566370614
datapath = '/home/soumitra/Solar_Flare_Project_Machine_Learning/Final' #JSOC files path definition
os.chdir(datapath)

#===========================================

def main():
   files_bz= sorted(glob.glob('*.Br.fits'))
   files_by= sorted(glob.glob('*.Bt.fits'))
   files_bx= sorted(glob.glob('*.Bp.fits'))
   files_bitmask= sorted(glob.glob('*.bitmap.fits'))
   files_mask=sorted(glob.glob('*.conf_disambig.fits'))
   files_los = sorted(glob.glob('*.magnetogram.fits'))
   ps=len(files_bz)
   print(ps)
   dim_nd_array= (ps, 24)
   ts_param= np.zeros(dim_nd_array)
   print(files_mask)
   print(files_bx)
   print(files_bitmask)
   print(files_los)
   print(files_bz)
#   ns=len(files_bz)-5
#   file = open("testfiless.txt", "w")

   for i in range(ps):
       file_bz=files_bz[i]
       print(file_bz)
       file_by=files_by[i]
       file_bx=files_bx[i]
       file_mask= files_mask[i]
       file_bitmask= files_bitmask[i]
       file_los=files_los[i]
       x = get_data(file_bz, file_by, file_bx, file_mask, file_bitmask, file_los)
       bz, by, bx, mask, bitmask, nx, ny, header, cdelt1, los = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]
       mean_vf, count_mask  = compute_abs_flux(bz, mask, bitmask, nx, ny, header, cdelt1)
       horiz      = compute_bh(bx, by, bz, mask, bitmask, nx, ny)
       bh = horiz[0]
       mean_gamma = compute_gamma(bx, by, bz, bh, mask, bitmask, nx, ny, header, cdelt1)
       total      = compute_bt(bx, by, bz, mask, bitmask, nx, ny)
       bt = total[0]
       lforc, lforcx, lforcy, lforcz, tforcx, tforcy, tforcz = compute_lorentz(bx, by, bz, mask, bitmask, nx, ny, header, cdelt1)
       mean_derivative_bt = computeBtderivative(bt, nx, ny, mask, bitmask)
       mean_derivative_bh = computeBhderivative(bh, nx, ny, mask, bitmask)
       mean_derivative_bz = computeBzderivative(bz, nx, ny, mask, bitmask)
       current                =  computeJz(bx, by, mask, bitmask, nx, ny)
       jz, derx, dery = current[0], current[1], current[2]
       mean_jz, us_i = computeJzmoments(jz, derx, dery, mask, bitmask, nx, ny, header, cdelt1, munaught)
       mean_alpha = computeAlpha(jz, bz, mask, bitmask, nx, ny, header, cdelt1)
       mean_ih, total_us_ih, total_abs_ih = computeHelicity(jz, bz, mask, bitmask, nx, ny, header, cdelt1)
       totaljz = computeSumAbsPerPolarity(jz, bz, mask, bitmask, nx, ny, header, cdelt1, munaught)
       potential = greenpot(bz, nx, ny)
       bpx, bpy  = potential[0], potential[1]
       meanpot, totpot = computeFreeEnergy(bx, by, bpx, bpy, nx, ny, header, cdelt1, mask, bitmask)
       meanshear_angle, area_w_shear_gt_45 = computeShearAngle(bx, by, bz, bpx, bpy, nx, ny, mask, bitmask)
       Rparam = computeR(los, nx, ny, cdelt1)
       ts_param[i, :] = [total_us_ih, lforc, totpot, total_us_ih, total_abs_ih, totaljz, mean_vf, tforcz, meanpot, Rparam, lforcz, area_w_shear_gt_45, meanshear_angle, mean_gamma, mean_derivative_bt,
                          mean_derivative_bz, mean_derivative_bh, mean_ih, tforcy, mean_jz, mean_alpha, tforcx, lforcy, lforcx] # arranged in Bobra's order; count_mask is considered as area_acr
#       file.write('%.5e   %.1f   %.5f   %.5f\n' % (mean_vf, count_mask, mean_jz, us_i))

   val_df = pd.DataFrame(ts_param)
#    df = pd.concat([ts_df, val_df], axis = 1)
   df=val_df
   df.columns = ['TOTUSJH', 'TOTBSQ', 'TOTPOT', 'TOTUSJZ', 'ABSNJZH', 'SAVNCPP', 'USFLUX', 'TOTFZ', 'MEANPOT', 'R_VALUE', 'EPSZ', 'SHRGT45', 'MEANSHR', 'MEANGAM', 'MEANGBT', 'MEANGBZ', 'MEANGBH',
                  'MEANJZH', 'TOTFY', 'MEANJZD', 'MEANALP', 'TOTFX', 'EPSY', 'EPSX']
    
   df.to_csv('/home/soumitra/Solar_Flare_Project_Machine_Learning/Final/24_params_total_final.csv', sep='\t', index = False)


def get_data(file_bz, file_by, file_bx, file_mask, file_bitmask, file_los):

    """function: get_data

    This function reads the appropriate data and metadata.
    """
    
    try:
        bz_map = sunpy.map.Map(file_bz)
    except:
        print("Could not open the bz fits file")
        sys.exit(1)

    try:
        by_map = sunpy.map.Map(file_by)
    except:
        print("Could not open the by fits file")
        sys.exit(1)

    try:
        bx_map = sunpy.map.Map(file_bx)
    except:
        print("Could not open the bx fits file")
        sys.exit(1)
    
    try:
        mask_map = sunpy.map.Map(file_mask)
    except:
        print("Could not open the mask fits file")
        sys.exit(1)

    try:
        bitmask_map = sunpy.map.Map(file_bitmask)
    except:
        print("Could not open the bitmask fits file")
        sys.exit(1)

    try:
        los_map = sunpy.map.Map(file_los)
    except:
        print("Could not open the magnetogram/los fits file")
        sys.exit(1)
  

    bx=bx_map.data
    by=by_map.data
    bz=bz_map.data
    mask=mask_map.data
    bitmask=bitmask_map.data
    los=los_map.data
    # get metadata
    header = bz_map.meta

    # convert cdelt1 from degrees to arcsec
    cdelt1 = (math.atan((header['rsun_ref']*header['cdelt1']*radsindeg)/(header['dsun_obs'])))*(1/radsindeg)*(3600.)

    # get dimensions
    nx     = bz.shape[1]
    ny     = bz.shape[0]

   # flip the sign of by
    by_flipped = -1.0*(np.array(by))
        
    return [bz, by_flipped, bx, mask, bitmask, nx, ny, header, cdelt1, los] 

#===========================================

def compute_abs_flux(bz, mask, bitmask, nx, ny, header, cdelt1):

    """function: compute_abs_flux

    This function computes the total unsigned flux in units of G/cm^2.
    It also returns the number of pixels used in this calculation in the keyword CMASK.
    
    To compute the unsigned flux, we simply calculate
       flux = surface integral [(vector Bz) dot (normal vector)],
            = surface integral [(magnitude Bz)*(magnitude normal)*(cos theta)].

    However, since the field is radial, we will assume cos theta = 1.
    Therefore, the pixels only need to be corrected for the projection.

    To convert G to G*cm^2, simply multiply by the number of square centimeters per pixel: 
       (Gauss/pix^2)(CDELT1)^2(RSUN_REF/RSUN_OBS)^2(100.cm/m)^2
       =Gauss*cm^2
    """

    count_mask = 0
    sum        = 0.0
        
    for j in range(ny):
        for i in range(nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(bz[j,i]):
                continue
            sum += abs(bz[j,i])
            count_mask += 1

    mean_vf     = sum*cdelt1*cdelt1*(header['rsun_ref']/header['rsun_obs'])*(header['rsun_ref']/header['rsun_obs'])*100.0*100.0
   
    return [mean_vf, count_mask]
       
#===========================================

def compute_bh(bx, by, bz, mask, bitmask, nx, ny):

    """function: compute_bh

    This function calculates B_h, the horizontal field, in units of Gauss.
    (The magnetic field has native units of Gauss since the filling factor = 1).
    """

    bh     = np.zeros([ny,nx])
    
    for j in range(ny):
        for i in range(nx):
            if (np.isnan(bx[j,i]) or np.isnan(by[j,i])):
                bh[j,i] = np.nan
                continue
            bh[j,i]     = np.sqrt(bx[j,i]*bx[j,i] + by[j,i]*by[j,i])
           
    return [bh]

#===========================================

def compute_lorentz(bx, by, bz, mask, bitmask, nx, ny, header, cdelt1):

    """function: compute_lorentz

    This function calculates lorentz force.
    (The magnetic field has native units of Gauss since the filling factor = 1).
    """

    bl     = np.zeros([ny,nx])
    bll     = np.zeros([ny,nx])
    bsx     = np.zeros([ny,nx])
    bsy     = np.zeros([ny,nx])
    sum=0.0
    suma=0.0
    sumb=0.0
    sumc=0.0
    
    for j in range(ny):
        for i in range(nx):
            if (np.isnan(bx[j,i]) or np.isnan(by[j,i]) or np.isnan(bz[j,i])):
                bl[j,i] = np.nan
                continue
            bl[j,i]     = bx[j,i]*bx[j,i] + by[j,i]*by[j,i] + bz[j,i]*bz[j,i]
            bll[j,i]     = bx[j,i]*bx[j,i] + by[j,i]*by[j,i] - bz[j,i]*bz[j,i]
            sum += bl[j,i]
            bsx[j,i]    = bx[j,i]*bz[j,i]
            bsy[j,i]    = by[j,i]*bz[j,i]
            suma += bsx[j,i]
            sumb += bsy[j,i]
            sumc += bll[j,i]

    lforc= sum
    lforcx= suma/sum
    lforcy= - sumb/sum
    lforcz = sumc/sum
    tforcx= - suma*cdelt1*cdelt1*(header['rsun_ref']/header['rsun_obs'])*(header['rsun_ref']/header['rsun_obs'])*100.0*100.0
    tforcy= sumb*cdelt1*cdelt1*(header['rsun_ref']/header['rsun_obs'])*(header['rsun_ref']/header['rsun_obs'])*100.0*100.0
    tforcz= sumc*cdelt1*cdelt1*(header['rsun_ref']/header['rsun_obs'])*(header['rsun_ref']/header['rsun_obs'])*100.0*100.0
       
    return [lforc, lforcx, lforcy, lforcz, tforcx, tforcy, tforcz]

#===========================================

def compute_gamma(bx, by, bz, bh, mask, bitmask, nx, ny, header, cdelt1):

    """function: compute_gamma

    This function computes the inclination of the horizontal field (relative to the radial field).

    """

    count_mask = 0
    sum        = 0.0
    
    for j in range(ny):
        for i in range(nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if ( np.isnan(bz[j,i]) or np.isnan(bh[j,i]) or bz[j,i] == 0 ):
                continue
            if ( bh[j,i] < 100 ):
                continue            
            sum += abs(math.atan(bh[j,i]/abs(bz[j,i])))*(180./np.pi)            
            count_mask += 1
                
    mean_gamma     = sum/count_mask
       
    return mean_gamma
    
#===========================================

def compute_bt(bx, by, bz, mask, bitmask, nx, ny):

    """function: compute_bt

    This function calculates B_t, the total field, in units of Gauss.
    (The magnetic field has native units of Gauss since the filling factor = 1).
    """

    bt     = np.zeros([ny,nx])    
    
    for j in range(ny):
        for i in range(nx):
            if (np.isnan(bx[j,i]) or np.isnan(by[j,i]) or np.isnan(bz[j,i])):
                bt[j,i] = np.nan                
                continue
            bt[j,i]     = np.sqrt(bx[j,i]*bx[j,i] + by[j,i]*by[j,i] + bz[j,i]*bz[j,i])
            
    return [bt]

#===========================================

def computeBtderivative(bt, nx, ny, mask, bitmask):

    """function: computeBtderivative

    This function computes the derivative of the total field.
    """

    count_mask = 0
    sum        = 0.0
    
    derx_bt   = np.zeros([ny,nx])
    dery_bt   = np.zeros([ny,nx])
    
    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1,nx-1):
        for j in range(0,ny):
           derx_bt[j,i]   = (bt[j,i+1] - bt[j,i-1])*0.5
               
    #brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0,nx):
        for j in range(1,ny-1):
           dery_bt[j,i]   = (bt[j+1,i] - bt[j-1,i])*0.5
    
    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero. 
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the mask and bitmask arrays.

    i=0
    for j in range(ny):
        derx_bt[j,i] = ( (-3*bt[j,i]) + (4*bt[j,i+1]) - (bt[j,i+2]) )*0.5
        
    i=nx-1
    for j in range(ny):
        derx_bt[j,i] = ( (3*bt[j,i]) + (-4*bt[j,i-1]) - (-bt[j,i-2]) )*0.5
    
    j=0
    for i in range(nx):
        dery_bt[j,i] = ( (-3*bt[j,i]) + (4*bt[j+1,i]) - (bt[(j+2),i]) )*0.5
    
    j=ny-1
    for i in range(nx):
        dery_bt[j,i] = ( (3*bt[j,i]) + (-4*bt[j-1,i]) - (-bt[j-2,i]) )*0.5

    # Calculate the sum only
    for j in range(1,ny-1):
        for i in range (1,nx-1):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if ( (derx_bt[j,i] + dery_bt[j,i]) == 0):
                continue
            if np.isnan(bt[j,i]):
                continue
            if np.isnan(bt[j+1,i]):
                continue
            if np.isnan(bt[j-1,i]):
                continue
            if np.isnan(bt[j,i-1]):
                continue
            if np.isnan(bt[j,i+1]):
                continue
            if np.isnan(derx_bt[j,i]):
                continue
            if np.isnan(dery_bt[j,i]):
                continue
            sum += np.sqrt( derx_bt[j,i]*derx_bt[j,i]  + dery_bt[j,i]*dery_bt[j,i]  )
            count_mask += 1

    mean_derivative_bt     = (sum)/(count_mask)

    return mean_derivative_bt

#===========================================

def computeBhderivative(bh, nx, ny, mask, bitmask):

    """function: computeBhderivative

    This function computes the derivative of the horizontal field.
    """

    count_mask = 0
    sum        = 0.0

    derx_bh   = np.zeros([ny,nx])
    dery_bh   = np.zeros([ny,nx])    
    
    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1,nx-1):
        for j in range(0,ny):
           derx_bh[j,i]   = (bh[j,i+1] - bh[j,i-1])*0.5
              
    #brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0,nx):
        for j in range(1,ny-1):
           dery_bh[j,i]   = (bh[j+1,i] - bh[j-1,i])*0.5
               
    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero. 
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the mask and bitmask arrays.

    i=0
    for j in range(ny):
        derx_bh[j,i] = ( (-3*bh[j,i]) + (4*bh[j,i+1]) - (bh[j,i+2]) )*0.5
        
    i=nx-1
    for j in range(ny):
        derx_bh[j,i] = ( (3*bh[j,i]) + (-4*bh[j,i-1]) - (-bh[j,i-2]) )*0.5
    
    j=0
    for i in range(nx):
        dery_bh[j,i] = ( (-3*bh[j,i]) + (4*bh[j+1,i]) - (bh[(j+2),i]) )*0.5
    
    j=ny-1
    for i in range(nx):
        dery_bh[j,i] = ( (3*bh[j,i]) + (-4*bh[j-1,i]) - (-bh[j-2,i]) )*0.5

    # Calculate the sum only
    for j in range(1,ny-1):
        for i in range (1,nx-1):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if ( (derx_bh[j,i] + dery_bh[j,i]) == 0):
                continue
            if np.isnan(bh[j,i]):
                continue
            if np.isnan(bh[j+1,i]):
                continue
            if np.isnan(bh[j-1,i]):
                continue
            if np.isnan(bh[j,i-1]):
                continue
            if np.isnan(bh[j,i+1]):
                continue
            if np.isnan(derx_bh[j,i]):
                continue
            if np.isnan(dery_bh[j,i]):
                continue
            sum += np.sqrt( derx_bh[j,i]*derx_bh[j,i]  + dery_bh[j,i]*dery_bh[j,i]  )
            count_mask += 1

    mean_derivative_bh     = (sum)/(count_mask)
    
    return mean_derivative_bh

#===========================================

def computeBzderivative(bz, nx, ny, mask, bitmask):

    """function: computeBzderivative

    This function computes the derivative of the vertical field.
    """

    count_mask = 0
    sum        = 0.0
    
    derx_bz   = np.zeros([ny,nx])
    dery_bz   = np.zeros([ny,nx])
        
    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1,nx-1):
        for j in range(0,ny):
           derx_bz[j,i]   = (bz[j,i+1] - bz[j,i-1])*0.5
              
    #brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0,nx):
        for j in range(1,ny-1):
           dery_bz[j,i]   = (bz[j+1,i] - bz[j-1,i])*0.5
              
    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero. 
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the mask and bitmask arrays.

    i=0
    for j in range(ny):
        derx_bz[j,i] = ( (-3*bz[j,i]) + (4*bz[j,i+1]) - (bz[j,i+2]) )*0.5
        
    i=nx-1
    for j in range(ny):
        derx_bz[j,i] = ( (3*bz[j,i]) + (-4*bz[j,i-1]) - (-bz[j,i-2]) )*0.5
    
    j=0
    for i in range(nx):
        dery_bz[j,i] = ( (-3*bz[j,i]) + (4*bz[j+1,i]) - (bz[(j+2),i]) )*0.5
    
    j=ny-1
    for i in range(nx):
        dery_bz[j,i] = ( (3*bz[j,i]) + (-4*bz[j-1,i]) - (-bz[j-2,i]) )*0.5

    # Calculate the sum only
    for j in range(1,ny-1):
        for i in range (1,nx-1):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if ( (derx_bz[j,i] + dery_bz[j,i]) == 0):
                continue
            if np.isnan(bz[j,i]):
                continue
            if np.isnan(bz[j+1,i]):
                continue
            if np.isnan(bz[j-1,i]):
                continue
            if np.isnan(bz[j,i-1]):
                continue
            if np.isnan(bz[j,i+1]):
                continue            
            if np.isnan(derx_bz[j,i]):
                continue
            if np.isnan(dery_bz[j,i]):
                continue
            sum += np.sqrt( derx_bz[j,i]*derx_bz[j,i]  + dery_bz[j,i]*dery_bz[j,i]  )
            count_mask += 1

    mean_derivative_bz     = (sum)/(count_mask)
    
    return mean_derivative_bz

#===========================================

def computeJz(bx, by, mask, bitmask, nx, ny):

    """function: computeJz

    This function computes the z-component of the current.

    In discretized space like data pixels, the current (or curl of B) is calculated as the integration
    of the field Bx and By along the circumference of the data pixel divided by the area of the pixel.

    One form of differencing the curl is expressed as:
    (dx * (Bx(i,j-1)+Bx(i,j)) / 2
    +dy * (By(i+1,j)+By(i,j)) / 2
    -dx * (Bx(i,j+1)+Bx(i,j)) / 2
    -dy * (By(i-1,j)+By(i,j)) / 2) / (dx * dy)

    To change units from Gauss/pixel to mA/m^2 (the units for Jz in Leka and Barnes, 2003),
    one must perform the following unit conversions:
    (Gauss)(1/arcsec)(arcsec/meter)(Newton/Gauss*Ampere*meter)(Ampere^2/Newton)(milliAmpere/Ampere), or
    (Gauss)(1/CDELT1)(RSUN_OBS/RSUN_REF)(1 T / 10^4 Gauss)(1 / 4*PI*10^-7)( 10^3 milliAmpere/Ampere), or
    (Gauss)(1/CDELT1)(RSUN_OBS/RSUN_REF)(0.00010)(1/MUNAUGHT)(1000.),
    where a Tesla is represented as a Newton/Ampere*meter.

    The units of total unsigned vertical current (us_i) are simply in A. In this case, we would have the following:
    (Gauss/pix)(1/CDELT1)(RSUN_OBS/RSUN_REF)(0.00010)(1/MUNAUGHT)(CDELT1)(CDELT1)(RSUN_REF/RSUN_OBS)(RSUN_REF/RSUN_OBS)
    = (Gauss/pix)(0.00010)(1/MUNAUGHT)(CDELT1)(RSUN_REF/RSUN_OBS)
    """

    count_mask = 0
    sum        = 0.0
    
    derx      = np.zeros([ny,nx])
    dery      = np.zeros([ny,nx])
    jz        = np.zeros([ny,nx])
    
    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1,nx-1):
        for j in range(0,ny):
           derx[j,i]      = (by[j,i+1] - by[j,i-1])*0.5
               
    #brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0,nx):
        for j in range(1,ny-1):
           dery[j,i]      = (bx[j+1,i] - bx[j-1,i])*0.5
           
    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero. 
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the mask and bitmask arrays.

    i=0
    for j in range(ny):
        derx[j,i] = ( (-3*by[j,i]) + (4*by[j,i+1]) - (by[j,i+2]) )*0.5
        
    i=nx-1
    for j in range(ny):
        derx[j,i] = ( (3*by[j,i]) + (-4*by[j,i-1]) - (-by[j,i-2]) )*0.5
    
    j=0
    for i in range(nx):
        dery[j,i] = ( (-3*bx[j,i]) + (4*bx[j+1,i]) - (bx[(j+2),i]) )*0.5
    
    j=ny-1
    for i in range(nx):
        dery[j,i] = ( (3*bx[j,i]) + (-4*bx[j-1,i]) - (-bx[j-2,i]) )*0.5

    # Calculate the sum only
    for j in range(1,ny-1):
        for i in range (1,nx-1):
            jz[j,i]     = (derx[j,i] - dery[j,i])
            
    return [jz, derx, dery]

#===========================================

def computeJzmoments(jz, derx, dery, mask, bitmask, nx, ny, header, cdelt1, munaught):

    """function: computeJzmoments

    This function computes moments of the vertical current.
    The mean vertical current density is in units of mA/m^2.
    The total unsigned vertical current is in units of Amperes.
    """

    count_mask = 0
    curl       = 0.0    
    us_i       = 0.0

    # Calculate the sum only
    for j in range(ny):
        for i in range(nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(jz[j,i]):
                continue
            if np.isnan(derx[j,i]):
                continue
            if np.isnan(dery[j,i]):
                continue        
            curl += (jz[j,i])*(1/cdelt1)*(header['rsun_obs']/header['rsun_ref'])*(0.00010)*(1/munaught )*(1000.)
            us_i += abs(jz[j,i])*(cdelt1/1)*(header['rsun_ref']/header['rsun_obs'])*(0.00010)*(1/munaught)
            count_mask += 1
            
    mean_jz     = curl/(count_mask)
        
    us_i        = (us_i)
    
    return [mean_jz, us_i]

#===========================================

def computeAlpha(jz, bz, mask, bitmask, nx, ny, header, cdelt1):

    """function: computeAlpha

    This function computes the twist parameter.

    The twist parameter, alpha, is defined as alpha = Jz/Bz. In this case, the calculation for alpha is weighted by Bz:

    numerator   = sum of all Jz*Bz
    denominator = sum of Bz*Bz
    alpha       = numerator/denominator

    The units of alpha are in 1/Mm
    The units of Jz are in Gauss/pix; the units of Bz are in Gauss.

    Therefore, the units of Jz/Bz = (Gauss/pix)(1/Gauss)(pix/arcsec)(arsec/meter)(meter/Mm), or
    = (Gauss/pix)(1/Gauss)(1/CDELT1)(RSUN_OBS/RSUN_REF)(10^6)
    = 1/Mm
    """

    alpha_total         = 0.0
    C                   = ((1/cdelt1)*(header['rsun_obs']/header['rsun_ref'])*(1000000.))
    total               = 0.0
    A                   = 0.0
    B                   = 0.0

    for j in range(ny):
        for i in range(nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(jz[j,i]):
                continue
            if np.isnan(bz[j,i]):
                continue
            if (jz[j,i] == 0):
                continue        
            if (bz[j,i] == 0):
                continue        
            A += jz[j,i]*bz[j,i]
            B += bz[j,i]*bz[j,i]
   
    #Determine the absolute value of alpha. The units for alpha are 1/Mm
    alpha_total         = ((A/B)*C)
    mean_alpha          = alpha_total
    
    return mean_alpha
   
#===========================================

def computeHelicity(jz, bz, mask, bitmask, nx, ny, header, cdelt1):

    """function: computeHelicity

    This function computes a proxy for the current helicity and various moments. 

    The current helicity is defined as Bz*Jz and the units are G^2 / m
    The units of Jz are in G/pix; the units of Bz are in G.
    Therefore, the units of Bz*Jz = (Gauss)*(Gauss/pix) = (Gauss^2/pix)(pix/arcsec)(arcsec/meter)
    = (Gauss^2/pix)(1/CDELT1)(RSUN_OBS/RSUN_REF)
    =  G^2 / m.
    """

    count_mask      = 0.0
    sum             = 0.0
    sum2            = 0.0
    
    for j in range(ny):
        for i in range (nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(jz[j,i]):
                continue
            if np.isnan(bz[j,i]):
                continue
            if (jz[j,i] == 0):
                continue        
            if (bz[j,i] == 0):
                continue        
            
            sum        += (jz[j,i]*bz[j,i])*(1/cdelt1)*(header['rsun_obs']/header['rsun_ref'])    #contributes to MEANJZH and ABSNJZH
            sum2       += abs(jz[j,i]*bz[j,i])*(1/cdelt1)*(header['rsun_obs']/header['rsun_ref']) # contributes to TOTUSJH
            count_mask += 1

    mean_ih          = sum/count_mask                                                               # Units are G^2 / m ; keyword is MEANJZH
    total_us_ih      = sum2                                                                         # Units are G^2 / m ; keyword is TOTUSJH
    total_abs_ih     = abs(sum)                                                                     # Units are G^2 / m ; keyword is ABSNJZH
         
    return [mean_ih, total_us_ih, total_abs_ih]

#===========================================

def computeSumAbsPerPolarity(jz, bz, mask, bitmask, nx, ny, header, cdelt1, munaught):

    """function: computeSumAbsPerPolarity

    This function computes the sum of the absolute value of the current per polarity. It is defined as follows:

    The sum of the absolute value per polarity is defined as the following:
    abs(sum(jz gt 0)) + abs(sum(jz lt 0)) and the units are in Amperes per arcsecond.
    The units of jz are in G/pix. In this case, we would have the following:
    Jz = (Gauss/pix)(1/CDELT1)(0.00010)(1/MUNAUGHT)(RSUN_REF/RSUN_OBS)(RSUN_REF/RSUN_OBS)(RSUN_OBS/RSUN_REF),
       = (Gauss/pix)(1/CDELT1)(0.00010)(1/MUNAUGHT)(RSUN_REF/RSUN_OBS)

    The error in this quantity is the same as the error in the mean vertical current.
    """
    
    count_mask      = 0.0
    sum1            = 0.0
    sum2            = 0.0

    for j in range(ny):
        for i in range (nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(jz[j,i]):
                continue
            if np.isnan(bz[j,i]):
                continue
            if (bz[j,i] > 0):
                sum1 += ( jz[j,i])*(1/cdelt1)*(0.00010)*(1/munaught)*(header['rsun_ref']/header['rsun_obs'])
            if (bz[j,i] <= 0):
                sum2 += ( jz[j,i])*(1/cdelt1)*(0.00010)*(1/munaught)*(header['rsun_ref']/header['rsun_obs'])
            count_mask += 1

    totaljz     = abs(sum1) + abs(sum2)
    
    return totaljz

#===========================================

def computeFreeEnergy(bx, by, bpx, bpy, nx, ny, header, cdelt1, mask, bitmask):
    """
    function: computeFreeEnergy

    This function computes the mean photospheric excess magnetic energy and total photospheric excess magnetic energy density.

    The units for magnetic energy density in cgs are ergs per cubic centimeter. The formula B^2/8*PI integrated over all space, dV
    automatically yields erg per cubic centimeter for an input B in Gauss. Note that the 8*PI can come out of the integral; thus,
    the integral is over B^2 dV and the 8*PI is divided at the end.

    Total magnetic energy is the magnetic energy density times dA, or the area, and the units are thus ergs/cm. To convert
    ergs per centimeter cubed to ergs per centimeter, simply multiply by the area per pixel in cm:
    erg/cm^3*(CDELT1^2)*(RSUN_REF/RSUN_OBS ^2)*(100.^2)
    = erg/cm(1/pix^2)
    """
    count_mask      = 0.0
    sum             = 0.0
    sum1            = 0.0
    
    for j in range(ny):
        for i in range (nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(bx[j,i]):
                continue
            if np.isnan(by[j,i]):
                continue
            sum  += ( ((bx[j,i] - bpx[j,i])*(bx[j,i] - bpx[j,i])) + ((by[j,i] - bpy[j,i])*(by[j,i] - bpy[j,i])) )*(cdelt1*cdelt1*(header['rsun_ref']/header['rsun_obs'])*(header['rsun_ref']/header['rsun_obs'])*100.0*100.0)
            sum1 += (  ((bx[j,i] - bpx[j,i])*(bx[j,i] - bpx[j,i])) + ((by[j,i] - bpy[j,i])*(by[j,i] - bpy[j,i])) )
            count_mask += 1

    # Units of meanpotptr are ergs per centimeter
    meanpot      = (sum1) / (count_mask*8.*np.pi)
        
    # Units of sum are ergs/cm^3, units of factor are cm^2/pix^2; therefore, units of totpotptr are ergs per centimeter
    totpot       = (sum)/(8.*np.pi)
        
    return [meanpot, totpot]

#===========================================

def computeShearAngle(bx, by, bz, bpx, bpy, nx, ny, mask, bitmask):
    """
    function: computeShearAngle

    This function computes the shear angle, or the angle between the potential field vector and the observed field vector, in degrees.
    """
    
    count_mask          = 0.0
    count               = 0.0
    dotproduct          = 0.0
    magnitude_potential = 0.0
    magnitude_vector    = 0.0
    sumsum              = 0.0
    shear_angle         = 0.0
 
    for j in range(ny):
        for i in range (nx):
            if ( mask[j,i] < 70 or bitmask[j,i] < 30 ):
                continue
            if np.isnan(bx[j,i]):
                continue
            if np.isnan(by[j,i]):
                continue
            if np.isnan(bz[j,i]):
                continue
            if np.isnan(bpx[j,i]):
                continue
            if np.isnan(bpy[j,i]):
                continue
          # for the values
            dotproduct            = (bpx[j,i])*(bx[j,i]) + (bpy[j,i])*(by[j,i]) + (bz[j,i])*(bz[j,i])
            magnitude_potential   = np.sqrt( (bpx[j,i]*bpx[j,i]) + (bpy[j,i]*bpy[j,i]) + (bz[j,i]*bz[j,i]))
            magnitude_vector      = np.sqrt( (bx[j,i]*bx[j,i])   + (by[j,i]*by[j,i])   + (bz[j,i]*bz[j,i]) )
            shear_angle           = math.acos(dotproduct/(magnitude_potential*magnitude_vector))*(180./np.pi)
            sumsum                += shear_angle
            count                 += 1
            if (shear_angle > 45):
                count_mask += 1

    # For mean 3D shear angle, area with shear greater than 45
    meanshear_angle     = (sumsum)/(count)
    
    # The area here is a fractional area -- the % of the total area. This has no error associated with it.
    area_w_shear_gt_45   = (count_mask/(count))*(100.0)

    return [meanshear_angle, area_w_shear_gt_45]

#===========================================

def computeR(los, nx, ny, cdelt1):
    """
    function: computeR

    This function computes the gradient-weighted neutral line length in Maxwells.
    """

    sum   = 0.0
    err   = 0.0
    sigma = 10.0/2.3548
    scale = int(round(2.0/cdelt1))

    # =============== [STEP 1] =============== 
    # bin the line-of-sight magnetogram down by a factor of scale using nearest-neighbor interpolation
    xvalues            = [range(0,nx,scale)]*int(round(ny/scale))
    # ideally the expression below should be (round(nx/scale)) for (len(xvalues[0])), but sometimes
    # this means that len(yvalues[0]) > or < len(xvalues[0]). thus it is better this way:
    yvalues            = [[i]*int(len(xvalues[0])) for i in range(0,ny,scale)]
    
    # add some conditions to make sure xvalues and yvalues are of the same dimension
    if len(xvalues) > len(yvalues):
        xvalues = xvalues[:int(len(yvalues))]
    if len(yvalues) > len(xvalues):
        yvalues = yvalues[:int(len(xvalues))]
    
    interp_coordinates = (xvalues,yvalues)    
    rim = scipy.ndimage.interpolation.map_coordinates(los,interp_coordinates,mode='nearest')

    # =============== [STEP 2] =============== 
    # identify positive and negative pixels greater than +/- 150 gauss
    # and label those pixels with a 1.0 in arrays p1p0 and p1n0

    nx1  = rim.shape[1]
    ny1  = rim.shape[0]
    p1p0 = np.zeros([ny1,nx1])
    p1n0 = np.zeros([ny1,nx1])

    for j in range(ny1):
        for i in range (nx1):
            if (rim[j,i] > 150):
                p1p0[j,i]=1.0
            else:
                p1p0[j,i]=0.0
            if (rim[j,i] < -150):
                p1n0[j,i]=1.0
            else:
                p1n0[j,i]=0.0

    # =============== [STEP 3] =============== 
    # smooth each of the negative and positive pixel bitmaps by convolving with a boxcar     

    # set up the convolution kernel
    boxcar_kernel = np.zeros([ny1,nx1])
    midpoint_ny1  = int(round(ny1/2))
    midpoint_nx1  = int(round(nx1/2))

    for j in range(midpoint_ny1,midpoint_ny1+3):
        for i in range(midpoint_nx1,midpoint_nx1+3):
            boxcar_kernel[j,i]=0.1111

    p1p = scipy.ndimage.filters.convolve(p1p0,boxcar_kernel)
    p1n = scipy.ndimage.filters.convolve(p1n0,boxcar_kernel)

    # =============== [STEP 4] =============== 
    # find the pixels for which p1p and p1n are both equal to 1. 
    # this defines the polarity inversion line

    p1 = np.zeros([ny1,nx1])
    for j in range(ny1):
        for i in range (nx1):
            if ((p1p[j,i] > 0.0) and (p1n[j,i] > 0.0)):
                p1[j,i]=1.0
            else:
                p1[j,i]=0.0
                
    # =============== [STEP 5] =============== 
    # convolve the polarity inversion line map with a gaussian
    # to identify the region near the plarity inversion line
    # the resultant array is called pmap

    pmap = scipy.ndimage.filters.gaussian_filter(p1,sigma,order=0)

    # =============== [STEP 6] =============== 
    # the R parameter is calculated

    for j in range(ny1):
        for i in range (nx1):
            if np.isnan(pmap[j,i]):
                continue
            if np.isnan(rim[j,i]):
                continue
            sum += pmap[j,i]*abs(rim[j,i])

    if (sum < 1.0):
        Rparam = 0.0
    else:
        Rparam = math.log10(sum)

    #return [Rparam]
    return Rparam
    
#=========================================== 

def greenpot(bz, nx, ny):
    """
    function: greenpot

    This function extrapolates the potential magnetic field using Green's functions.
    The underlying assuption of a potential field is that it is Maxwell-stress free.
    The monopole depth is 0.01 pixels.
    """
    print ('Calculating the potential field. This takes a minute.')

    nnx = nx
    nny = ny
    
    # define the monopole depth, dz
    dz = 0.001

    # malloc some arrays
    pfpot      = np.zeros([nny,nnx])
    rdist      = np.zeros([nny,nnx])
    bztmp      = np.zeros([nny,nnx])
    bxp        = np.zeros([nny,nnx])
    byp        = np.zeros([nny,nnx])

    # substitute zeros for nans in bz data
    for iny in range(nny):
        for inx in range(nnx):
            if np.isnan(bz[iny,inx]):
                bztmp[iny,inx] = 0.0
            else:
                bztmp[iny,inx] = bz[iny,inx]

    rdd  = 0.0
    rdd1 = 0.0
    rdd2 = 0.0
    for iny in range(nny):
        for inx in range(nnx):
            rdd1  = float(inx)
            rdd2  = float(iny)
            rdd   = rdd1 * rdd1 + rdd2 * rdd2 + dz * dz
            rdist[iny,inx] = 1.0/(np.sqrt(rdd))

    iwindow = 0
    if (nnx > nny):
        iwindow = nnx
    else:
        iwindow = nny

    rwindow = float(iwindow)
    rwindow = rwindow * rwindow + 0.01 # must be square    
    rwindow = 1.0e2                    # limit the window size to be 10.    
    rwindow = np.sqrt(rwindow)
    iwindow = int(rwindow)

    for iny in range(nny):
        for inx in range(nnx):
            if np.isnan(bz[iny,inx]):
                pfpot[iny,inx] = 0.0
            else:
                sum = 0.0
                j2s = iny - iwindow
                j2e = iny + iwindow
                if (j2s < 0):
                    j2s = 0
                if (j2e > nny):
                    j2e = nny
                i2s = inx - iwindow
                i2e = inx + iwindow
                if (i2s < 0):
                    i2s = 0
                if (i2e > nnx):
                    i2e = nnx
                for j2 in range(j2s,j2e):
                    for i2 in range(i2s,i2e):
                        val1 = bztmp[j2,i2]
                        di = abs(i2 - inx)
                        dj = abs(j2 - iny)
                        sum = sum + val1 * rdist[dj,di] * dz
                pfpot[iny,inx] = sum
    
    for iny in range(1,nny-1):
        for inx in range(1,nnx-1):
            bxp[iny,inx] = -(pfpot[iny,inx+1] - pfpot[iny,inx-1])*0.5
            byp[iny,inx] = -(pfpot[iny+1,inx] - pfpot[iny-1,inx])*0.5

    return [bxp, byp]

#===========================================
  
if __name__ == "__main__":
    main()
    
__author__ = 'Soumitra Hazra'
