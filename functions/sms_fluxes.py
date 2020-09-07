##### Repository for the calculations, mostly related to submesoscale heat fluxes
##### Isabelle Giddy 2020

###########################################################

# All data should be gridded to equal horizontal distances to standardise gradients in density

# Python Modules

import numpy as np
import matplotlib.pyplot as plt

# Some standard constants are required

#g=9.8   --- Acceleration due to gravity
#cp=4000   --- Heat capacity of seawater

#################################################################################################

# Coriolis force

def coriolis_acc(lat,omega=0.729*10**-4):
    return 2*omega*np.sin(lat/360.*2*np.pi)

# Mixed Layer Depth

def calc_mld(density,depth, ref_depth=10,threshold=0.03):
    """
    Calculated Mixed Layer Depth
    Default threshold is 0.03 referenced to 10m
    Data should be in format depth,distance
    """

    mld=np.ndarray(len(density[1,:]))
    for i in range(len(density[1,:])):
        try: mld[i]=(depth[(np.abs((density[:,i]-density[ref_depth,i ]))>=threshold)].min())
        except ValueError:  #raised if `y` is empty.
            mld[i]=(np.nan)
    return mld



# Alternative mld based on n2
def cal_mldn2(density,depth,n2,ref_depth=10,threshold=0.05):
    mld=np.ndarray(len(density[1,:]))
    for i in range(len(density[1,:])):
        try:
            drange = (depth[(np.abs((density[:,i]-density[ref_depth,i ]))>=threshold)].min())
            mld[i]=depth[depth<=drange][n2[depth<=drange,i]==np.nanmax(n2[depth<=drange,i])]#

        except ValueError:  #raised if `y` is empty.
            mld[i]=(np.nan)
    return mld
#################################################################################################

# thermal expansion co-efficient, saline contraction coefficient

def alphabeta(SA,CT,depth):
    import gsw as gsw

    _,y=np.meshgrid(SA[1,:],depth)
    _, alpha, beta = gsw.rho_alpha_beta(SA,CT,y) 
    return alpha,beta

#################################################################################################

# Buoyancy and buoyancy gradient

def calc_buoyancy(density,SA,CT,alpha,beta,mld,dx=1000,po=1027,g=9.8):
    """
    Calculates buoyancy, buoyancy gradient, mean buoyancy gradient in the mixed layer
    """

    by = g * (1-density / po)   # Buoyancy

    by_S = g * (1-beta*SA/po)

    by_T = g * (1-alpha*CT/po)


    #Raw buoyancy gradient
    bgrad = (np.diff(by,axis=1)/dx)

    #Buoyancy gradient in the middle of the mixed layer
    bxml=np.ndarray(len(density[1,:-1]))
    for i in range(len(density[1,:-1])):
        bxml[i]=(np.nanmean(bgrad[:np.int8(mld[i])-15,i],0))

    return by,by_S,by_T,bgrad,bxml

#################################################################################################

# Glider Trajectory

def calculate_initial_compass_bearing(pointA, pointB):
    """
    Calculates the bearing between two points.
    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees
    :Returns Type:
      float
    """
    import math
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
            * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(y, x)
   # initial_bearing = np.atan2(x, y)


    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    #compass_bearing = (initial_bearing + 360) % 360

    return initial_bearing

#def glider_bearing(lon,lat):
#    import math
#    lon=np.deg2rad(lon)
##    lat=np.deg2rad(lat)
#    bearing = np.ndarray(len(lon[:-1]))
    #bearing=[]#
#    for i in range(len(lon[:-1])):
#        dlon=np.deg2rad(lon[i+1]-lon[i])

#        deltaX = math.cos(lat[i])*math.sin(lat[i+1])-math.sin(lat[i])*math.cos(lat[i+1])*math.cos(dlon)
#        deltaY = math.sin(dlon) * math.cos(lat[i+1])
        #convert to degrees
#        bearing[i]=(math.atan2(deltaX, deltaY))* (180/math.pi) # Compute such that 0 degrees is east
#    return bearing
    #normalize to compass headings
    #bearing = (bearing + 180) % 360

def calc_glider_traj(lat,lon):
    bearing=[]
    for i in range(len(lon[:-1])):
        bearing.append(calculate_initial_compass_bearing((lat[i],lon[i]),(lat[i+1],lon[i+1])))
    #bearing=np.array(bearing)
    return bearing

#################################################################################################


# Function to compute the difference between two angles

def calculateDifferenceBetweenAngles(bearing1, bearing0):
    norm1 = (bearing1 + 360) % 360
    norm2 = (bearing0 + 360) % 360
    
  
  #  norm1 = bearing1
  #  norm2=bearing0

    if (norm1)<(norm2):

        diffangle = (norm1 - norm2) + 180
        diffangle = (diffangle / 180.0)
        diffangle = ((diffangle - np.floor( diffangle )) * 360.0) - 180
    else:
        diffangle = (norm2 - norm1) + 180
        diffangle = (diffangle / 180.0)
        diffangle = ((diffangle - np.floor( diffangle )) * 360.0) - 180
 
    return diffangle


#################################################################################################

# Compute wind stress and direction

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def calc_wind(u10_wind,v10_wind):

    import airsea
    from metpy.calc import wind_direction
    tau_x=airsea.windstress.stress(u10_wind,z=10,drag='largepond',rho_air=1.22,Ta=10.)
    tau_y=airsea.windstress.stress(v10_wind,z=10,drag='largepond',rho_air=1.22,Ta=10.)


    tx = np.where(u10_wind>0,tau_x,tau_x*-1)  # Add directional signs back to wind stress vectors
    ty = np.where(v10_wind>0,tau_y,tau_y*-1)

    tau = np.sqrt(tx**2+ty**2)

    [tau2,theta]=cart2pol(u10_wind,v10_wind)   #to get winds oriented 0 for Easterly. 
    wind_dir=np.rad2deg(theta)

    return tau, tx, ty, wind_dir

#############################################################
# Plot histogram of wind orientation

def plot_wind_orientation(wind_dir,plot_title):

    normalise = (wind_dir+ 360) % 360   # First normalise winds from -180-180 to 0-360

    ax = plt.subplot(111, projection='polar')
    histogram, bins = np.histogram(np.deg2rad(normalise), bins=25)
    bin_centers = 0.5*(bins[1:] + bins[:-1])
    ax.bar(bin_centers, histogram,label="Wind Orientation",color='lightblue',bottom=0.0,alpha=0.8,edgecolor='tab:blue')
    ax.grid(alpha=0.2) 
    ax.yaxis.get_major_locator().base.set_params(nbins=5)

    ax.set_xlabel('{}'.format(plot_title))

    return ax

#################################################################################################

# Interpolate winds to glider time step

def glider_compat(var_time,glider_time,var,time_raw=False):

    if time_raw == False:
        #convert time
        var_time = np.int64(var_time)//10**9 * 10**9
        sg_time = np.int64(glider_time)//10**9 * 10**9
            #interp era5 time to sg time
        var_interp = np.interp(sg_time,var_time,var.squeeze())

    else:
        sg_time = np.int64(glider_time)//10**9 * 10**9

        #interp era5 time to sg time
        var_interp = np.interp(sg_time,var_time,var.squeeze())

    return var_interp

#################################################################################################


# Rotate winds to glider trajectory

def rotate_winds(wind_dir,glider_dir,tau_x,tau_y):
    
    bearing1 = wind_dir
    bearing0 = glider_dir
    x = np.abs(tau_x)
    y = np.abs(tau_y)

    angle=np.ndarray(len(bearing1))
    for k in range(len(bearing1)):
        angle[k]=(calculateDifferenceBetweenAngles(bearing1[k],bearing0[k]))
    angle=180-angle
    theta=np.deg2rad(angle)
    

    Rtx= x*np.cos(theta)+(y*np.sin(theta))  
    return Rtx,angle

#################################################################################################


# EQUIVALENT HEAT FLUX OF MIXED LAYER EDDIES

def calc_qmle(buoyancy_gradient_mld,mld,alpha,f=1e-4,cp=4000,po=1027,g=9.8):
    """
    Calculates Qmle based on Fox-Kemper 2008 
    Strong lateral gradients provide a resevoir of potential energy which can be released b ageostrophic overturning circulations as a result of 
    ageostrophic baroclinic instabilities.
    The restratificatoin by ABI is slower than by Symmetric Instabilities but faster than mesoscale variations
    Here, ABI is computed as an equivalent heat flux 
    """
    qmle=0.06*((buoyancy_gradient_mld**2*(mld)**2)/np.abs(f))*((cp*po)/(alpha*g))
    return qmle

#################################################################################################


# EKMAN BUOYANCY FLUX

def calc_ebf(buoyancy_gradient_mld,wind_dir,glider_dir,tau_x,tau_y,f,alpha,cp=4000,g=9.8):

    """
    Calculates the wind force on the ocean over fronts
    Downfront windstress, or wind stress in the direction of geostrophic shear produces an Ekman transport that advects less buoyant over
    more buoyant water. An upfront windstress act to restratify the mixed layer. 
    Theory from Thomas, 2005; Thompson et al., 2016

    """

    rotated_wind_component,angle=rotate_winds(wind_dir,glider_dir,tau_x,tau_y)
    ebf=(-(buoyancy_gradient_mld)*np.array((rotated_wind_component/f))*(cp/(alpha*g)))

    return ebf, rotated_wind_component,angle


##################################################################################################

def calc_n2(SA,CT,rho,depth,alpha,g=9.8,po=1027,gsw=True):
    """
    gsw implementation of n2, 
    decomposed into temp and salt components
    """
    if gsw==True:
        # From gsw package
        import gsw
        n2=gsw.Nsquared(SA,CT,depth,axis=0) 
        n2_t=(-g*alpha*np.diff(CT,axis=0))   # Temperature component
        n2_s = n2-n2_t # Salt component
    else:
        # From buoyancy
        by = g * (1-rho / po)   # Buoyancy

        n2 = np.diff(by*-1,axis=0)
        n2_t=(-g*alpha[:-1,:]*np.diff(CT,axis=0))   # Temperature component
        n2_s = n2-n2_t # Salt component

    return n2,n2_t,n2_s


##################################################################################################

def calc_Lr(rho,mld,f,g=9.8,po=1027.):
    """
    Calculates in internal rossby radius of deformation 
    based on Timmermans and Windsor (2013)
    Generally defined as Lr = NH/f 
    """
    n2ml=np.ndarray(len(rho[1,:-1]))
    for i in range(len(rho[1,:-1])):
        n2ml[i]=-(g/po)*((rho[15,i]-rho[np.int8(mld[i])+15,i])/mld[i])
    Lr=(np.sqrt(n2ml)*mld[:-1])/f

    return Lr

##################################################################################################



# Ertel Potential Vorticity

def calc_ertelPV(n2, bx, rel_vorticity, g=9.8,f=-1e-4):
    """
    As in Thomas et al., 2013; Thompson et al., 2016
    The assumption is made that flow is in geostrophic balance

    PV can be used as a diagnostic to indicate the suscepptibility of a flow to instability e.g., Hoskins, 1974; Thomas et al., 2008
    When PV takes the opposite sign of f, the flow is inherently unstable (Thomas et al., 2013)

    PV can be decomposed into two components:
    qvert: The vertical component of absolute vorticity and vertical stratificaton (N2)
    qbc: The horizontal relative vorticity (vertical shear) and M2, the horizontal buoyancy gradient 

    """

    # vertical component

    qvert = (f+rel_vorticity)*n2

    # baroclinic component
    qbc = -bx**2/f

    # Ertel PV

    ertelPV = qvert + qbc

    # If PV is unstable
    fq = ertelPV*f # fq > 0 stable

    return  ertelPV, qvert, qbc, fq

##################################################################################################

# Compute Richardson and Rossby Number

def calc_ri(n2,bx,f):
    """
   Calculates the Non-Dimensional modified Richardson Number
   Thermal wind balance is assumed
   Gives a measure of vertical shear 
   Ri characterises the dynamic regime and can be interpreted as the steepness of isopycnal slopes relatively to N/f.
   Ri >> 1: QG regime, steepness of isopycnal slope is small
   Ri ~ 1 (<4), Ageostrophic regime, steepness of isopycnal slope is large
   Also used to characterise instabilities as in Thomas et al., 2013
    """

    Ri = (f**2*n2)/np.abs(bx)**2
    modRi = np.arctan(-Ri**-1)
    return Ri, modRi

def calc_ro(f,vorticity):
    """
    The modified Rossby Number gives a measure of horizontal shear 
    """

    Ro = vorticity/f
    modRo = np.arctan(-1-Ro)

    return Ro, modRo







