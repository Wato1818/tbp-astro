import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import scipy.stats as st
np.seterr(divide = 'ignore')
from cycler import cycler
from pandas import *
from astropy.coordinates import SkyCoord
import pandas as pd
from astropy.table import Table
import astropy.units as u
import mpl_toolkits.mplot3d.axes3d as axes3d
import statistics 
import itertools as it
import sys
import os

# catalog_path=os.path.dirname(__file__)
# print(catalog_path+"\n\n\n\n")
catalog_orion = read_csv(os.path.dirname(__file__)+"/Catalogs/ONC_gaiaedr3&apogee_cone_0.5deg.csv")
parallax_min=0.5#PARAMETER 1
parallax_max=4  #PARAMMETER 2

catalog_orion = catalog_orion.loc[(catalog_orion['parallax']>=parallax_min)]
catalog_orion = catalog_orion.loc[catalog_orion['parallax']<=parallax_max]
catalog_orion['parallax']=catalog_orion['parallax'].apply(lambda x: round(x,2))
catalog_orion.reset_index(inplace=True)

#ORION CLUSTER DATA
ori_center = SkyCoord("5h35m15.68s -5d23m40s", frame='icrs', unit = "deg")
ori_dist = 403 # [pc]
par_ori = 2.48/1000 #Paralaje [arcsec]
ori_ra = 83.81533333*0.0174533# [rad] 
ori_dec = -5.39444444*0.0174533 # [rad]
v_rad_ori = 21.8 #Velocidad radial [km/s]
prop_ra_ori =  1.96 #Proper Motion [mas/yr]
prop_dec_ori = -0.77 #Proper Motion [mas/yr]
kappa = 4.74 #Convertion Factor mas/yr -> km/s (d>1kpc)


#ORION STARS DATA
ra_stars = catalog_orion["ra"]*0.0174533 #[rad]
ra_err_stars = catalog_orion["ra_error"]*0.0174533 #[rad]
dec_stars = catalog_orion["dec"]*0.0174533 #[rad] 
dec_err_stars = catalog_orion["dec_error"]*0.0174533 #[rad]

# pmra to pmra 2
prop_ra_stars = catalog_orion["pmra_2"]*np.cos(dec_stars) #[mas/yr]
prop_dec_stars = catalog_orion["pmdec"] #[mas/yr]
parallax_stars = catalog_orion["parallax"]/1000 #[arcsec]
parallax_err_stars = catalog_orion["parallax_error"]/1000 #[arcsec]
ruwe_orion = catalog_orion["ruwe"]
rad_vel = catalog_orion["dr2_radial_velocity"]
rad_vel_err = catalog_orion["dr2_radial_velocity_error"]
sum(rad_vel_err[rad_vel_err>0])/len(rad_vel_err[rad_vel_err>0])


#ADAPTING ORION STARS TO XY POSITIONS
x_ref_stars = (np.cos(dec_stars)*np.sin(ra_stars - ori_ra)) #~[rad]
y_ref_stars = (np.sin(dec_stars)*np.cos(ori_dec) - np.cos(dec_stars)*np.sin(ori_dec)*np.cos(ra_stars - ori_ra)) #~[rad]

#DISTANCE BETWEEN STARS
dist_xy = np.sqrt(x_ref_stars**2 + y_ref_stars**2) #[rad], orthographic distance
dist_xy_mas=dist_xy*10**3*3600/0.0174533 #~[mas] The proper motions are converted to radians really, this is usefull to tracebacktime.  

mu_ra_per_ori = (ra_stars-ori_ra)*(prop_dec_ori*np.sin(ori_dec) - (v_rad_ori*par_ori)/kappa*np.cos(ori_dec)) 
mu_dec_per_ori = -(ra_stars-ori_ra)*prop_ra_ori*np.sin(ori_dec) - (dec_stars-ori_dec)*(v_rad_ori*par_ori)/kappa

prop_obs_ra_ori = prop_ra_stars - prop_ra_ori
prop_obs_dec_ori = prop_dec_stars - prop_dec_ori
prop_rest_ra_ori = -(prop_obs_ra_ori - mu_ra_per_ori)
prop_rest_dec_ori = prop_obs_dec_ori - mu_dec_per_ori

prop_x_ori = prop_rest_ra_ori*np.cos(ra_stars - ori_ra) - prop_rest_dec_ori*np.sin(dec_stars)*np.sin(ra_stars - ori_ra) #[mas/yr]
prop_y_ori = prop_rest_ra_ori*np.sin(ori_dec)*np.sin(ra_stars - ori_ra) + prop_rest_dec_ori*(np.cos(dec_stars)*np.cos(ori_dec) 
+ np.sin(dec_stars)*np.sin(ori_dec)*np.cos(ra_stars - ori_ra)) #[mas/yr]

pm_xy_ori = np.sqrt(prop_x_ori**2 + prop_y_ori**2) #[mas/yr]
pm_xy_ori_kms=pm_xy_ori/(par_ori*1000)*4.74

interaction_orion = read_csv(os.path.dirname(__file__)+"/Interactions_catalog2.csv")
star_s = interaction_orion["star_s"]
star_b = interaction_orion["star_b"]
vel_mod_s = interaction_orion["vel_mod_s"]
vel_mod_b = interaction_orion["vel_mod_b"]
angle_sb = interaction_orion["angle"]
interaction_orion[(star_s==175) | (star_b==175)]
parallax_s = interaction_orion["parallax_s"]
parallax_b = interaction_orion["parallax_b"]
ruwe_s = interaction_orion["ruwe_s"]
ruwe_b = interaction_orion["ruwe_b"]


h=int(sys.argv[1])

n=(star_s[(star_s==h)|(star_b==h)].append(star_b[(star_s==h)|(star_b==h)])) #Star Number in the Catalog

fig, ax = plt.subplots()
m = h

for i in n:   
    x_pos = x_ref_stars[i] #Pos.x Runaways
    y_pos = y_ref_stars[i] #Pos.y Runaways
    x_direct = prop_x_ori[i] #PropMot Ra Runaways
    y_direct = prop_y_ori[i] #PropMot Dec Runaways
    x_pos_center = x_ref_stars[m] #Pos.x Main Binary
    y_pos_center = y_ref_stars[m] #Pos.y Main Binary
    x_direct_center = prop_x_ori[m] #PropMot Ra Main Binary
    y_direct_center = prop_y_ori[m] #PropMot Dec Main Binary
    
    #Plot Runaways Velocitie Data
    ax.quiver(x_pos, y_pos, x_direct, y_direct, width = 0.005, scale= 100, zorder = 1)
    #Plot the positions of the runaways in the Loop
    plt.scatter(x_ref_stars[i], y_ref_stars[i], marker = ".", s = 100)
    
    #Plot Main Binary Velocitie Data
    ax.quiver(x_pos_center, y_pos_center, x_direct_center, y_direct_center, width = 0.005, scale = 100)

    # plt.text(x_ref_stars[i], y_ref_stars[i], s = i, fontsize = 10) #Star number
    # plt.scatter(x_ref_stars[m], y_ref_stars[m], marker = ".", s = 100)
    plt.text(x_ref_stars[m], y_ref_stars[m], s = m, fontsize = 10)

ax.set_title('Runaways Velocities Orion')

#Plot the positions of Trapezium in the Loop
plt.scatter(0, 0, marker = "o",color ="green", s=10, zorder=2)

plt.xlabel("x")
plt.ylabel("y")
# plt.xlim(-0.01,0.01)
# plt.ylim(-0.01,0.01)
plt.savefig(os.path.dirname(__file__)+'/VelocityOrionStars2.jpg', dpi=300)
plt.grid()
# plt.show()
os.exit(0)