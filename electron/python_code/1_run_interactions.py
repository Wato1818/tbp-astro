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

#LOADING CATALOG
catalog_orion = read_csv("./Catalogs/ONC_gaiaedr3&apogee_cone_0.5deg.csv")
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


run_combinations = []
list = range(0,len(catalog_orion["source_id"]))
for c in it.combinations((list), 2):
    run_combinations.append(c)
# len(run_combinations)


star_s = []
star_b = []
cross_point_stars_x = []
cross_point_stars_y = []
TBL_S=[]#* u.ys
TBL_B=[]#* u.ys
angle_bet_stars = []
VelMod_S=[]#* u.m / u.s
VelMod_B=[]#* u.m / u.s
parallax_S = []
parallax_B = []
parallax_err_S = []
parallax_err_B = []
ruwe_S = []
ruwe_B = []

#MAGIC?
print("Doing the heavy stuff...")
def module(a,b):
    c = np.sqrt(a**2 + b**2)
    return c

for i in range(0,len(run_combinations)):
    s = run_combinations[i][0]
    b = run_combinations[i][1]
    
    #Position Single
    pos_s_x, pos_s_y = (x_ref_stars[s],y_ref_stars[s])
    #Velocity Single
    v_s_x, v_s_y = prop_x_ori[s], prop_y_ori[s]
    #Velocity Modules Single
    v_s_mod = module(v_s_x, v_s_y)
    #Norm Single
    nor_s_x, nor_s_y = (v_s_x/v_s_mod,v_s_y/v_s_mod)
    
    #Position Binary
    pos_b_x, pos_b_y =  (x_ref_stars[b],y_ref_stars[b])
    #Velocity Binary
    v_b_x, v_b_y = prop_x_ori[b], prop_y_ori[b]
    #Velocity Modules Binary
    v_b_mod = module(v_b_x,v_b_y)
    #Norm Binary
    nor_b_x, nor_b_y = (v_b_x/v_b_mod,v_b_y/v_b_mod)
    if v_s_mod>0 and v_b_mod>0:
        #Difference in positions of the Stars
        dif_pos_x = pos_s_x - pos_b_x
        dif_pos_y = pos_s_y - pos_b_y
        
        #Angle for the runaways
        angle_run = np.arccos((v_s_x*v_b_x + v_s_y*v_b_y)/(v_s_mod*v_b_mod))
        
        #Matrix of the equation system
        matrix = np.array([[nor_b_x,-nor_s_x],[nor_b_y, -nor_s_y]])
        matrix_res = np.array([dif_pos_x,dif_pos_y])
        matrix_solve = np.linalg.solve(matrix, matrix_res)
        
        if matrix_solve[0] < 0 and matrix_solve[1]<0:
            cross_point_sx = (matrix_solve[1]*nor_s_x + pos_s_x)
            cross_point_sy = (matrix_solve[1]*nor_s_y + pos_s_y)
            
            cross_point_bx = (matrix_solve[0]*nor_b_x + pos_b_x)
            cross_point_by = (matrix_solve[0]*nor_b_y + pos_b_y)
            
            traceback_line_b = module(pos_b_x 
            - cross_point_bx , pos_b_y - cross_point_by)/v_b_mod
            traceback_line_s = module(pos_s_x 
            - cross_point_sx , pos_s_y - cross_point_sy)/v_s_mod          
            
            if abs(traceback_line_b-traceback_line_s)<traceback_line_b/100*20 and angle_run>1.0:
                
                star_s.append(run_combinations[i][0])
                star_b.append(run_combinations[i][1])
                cross_point_stars_x.append(cross_point_sx)
                cross_point_stars_y.append(cross_point_sy)
                TBL_S.append(traceback_line_s*10**3*3600/0.074533) #[yr]
                TBL_B.append(traceback_line_b*10**3*3600/0.074533) #[yr]
                VelMod_S.append(pm_xy_ori_kms[s]) # [km/s]
                VelMod_B.append(pm_xy_ori_kms[b]) # [km/s]
                angle_bet_stars.append(angle_run) #[rad]
                parallax_S.append(parallax_stars[s]*10**3) #[mas]
                parallax_err_S.append(parallax_err_stars[s]*10**3) #[mas] 
                parallax_B.append(parallax_stars[b]*10**3) #[mas]
                parallax_err_B.append(parallax_err_stars[b]*10**3) #[mas]
                ruwe_S.append(ruwe_orion[s]) 
                ruwe_B.append(ruwe_orion[b])
                
T1 = Table([star_s,star_b, cross_point_stars_x, cross_point_stars_y, TBL_S, TBL_B, VelMod_S, VelMod_B, angle_bet_stars, parallax_S, parallax_err_S, parallax_B, parallax_err_B, ruwe_S, ruwe_B],names=("star_s","star_b", "cross_point_x", "cross_point_y", "TBL_S", "TBL_B", "vel_mod_s", "vel_mod_b","angle", "parallax_s","parallax_err_s", "parallax_b","parallax_err_b", "ruwe_s", "ruwe_b"),meta={'name': 'Table_three_body_interactions'})
T1.write('Interactions_catalog2.csv',format='ascii.csv', overwrite=True)

print("Finished")