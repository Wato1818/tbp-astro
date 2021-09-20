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

print(sys.argv[1])


interaction_orion = read_csv("Interactions_catalog2.csv")
star_s = interaction_orion["star_s"]
star_b = interaction_orion["star_b"]
vel_mod_s = interaction_orion["vel_mod_s"]
vel_mod_b = interaction_orion["vel_mod_b"]
angle_sb = interaction_orion["angle"]
# interaction_orion[(star_s==175) | (star_b==175)]
parallax_s = interaction_orion["parallax_s"]
parallax_b = interaction_orion["parallax_b"]
ruwe_s = interaction_orion["ruwe_s"]
ruwe_b = interaction_orion["ruwe_b"]

#SCORE CALCULATION
scores=[]
scores_id=[]

#---PARAMETERS
bs_max_velocity=999 #Binary Star maximum velocity
interactions_mode=True #if this value is True, the variable below works as a minimum. Else, the number of interactions must be equal to the value below.
bs_rule_interactions=10 #Binary Star rule interactions
#---PARAMETERS

for i in range(0,len(interaction_orion)):
    vel_mod_s=interaction_orion.iloc[i]['vel_mod_s']
    vel_mod_b=interaction_orion.iloc[i]['vel_mod_b']
    s=interaction_orion.iloc[i]['star_s']
    b=interaction_orion.iloc[i]['star_b']
    score = 0 #If the binary star does not accomplish the minimum requirements it gets a score of 0
    bs_id = s

    if((vel_mod_b<vel_mod_s) & (vel_mod_b<=bs_max_velocity)): #Filtering by binary star maximum velocity
        n_interactions=star_b[(star_b==b)|(star_s==b)].count()        
        if(interactions_mode):
            if(n_interactions>=bs_rule_interactions):
                vel_sum=sum(interaction_orion.vel_mod_s[(star_b==b)])
                score=vel_sum/n_interactions
                bs_id=b
        elif(n_interactions==bs_rule_interactions):
            vel_sum=sum(interaction_orion.vel_mod_s[(star_b==b)])
            score=vel_sum/n_interactions
            bs_id=b

    elif((vel_mod_s<vel_mod_b) & (vel_mod_s<=bs_max_velocity)):
        n_interactions=star_s[(star_s==s)|(star_b==s)].count()
        if(interactions_mode):
            if(n_interactions>=bs_rule_interactions):
                vel_sum=sum(interaction_orion.vel_mod_b[(star_s==s)])
                score=vel_sum/n_interactions
                bs_id=s
        elif(n_interactions==bs_rule_interactions):
            vel_sum=sum(interaction_orion.vel_mod_b[(star_s==s)])
            score=vel_sum/n_interactions
            bs_id=s
        

    scores.append(score)
    scores_id.append(bs_id)

test_interactions = interaction_orion
test_interactions["star_score"]=scores
test_interactions["star_score_ID"]=scores_id
test_interactions = test_interactions.sort_values(axis='index',kind='quicksort',by="star_score", ascending=False)

test_interactions.to_csv('scores.csv',encoding='utf-8')



#DISPLAY OPTIONS
# bs_min_score=-1
# bs_max_score=18 #Max 33
# max_to_display=5000
# allow_duplicate=False
# if(allow_duplicate):
#    filtered_interactions=test_interactions[(test_interactions["star_score"]>bs_min_score) & (test_interactions["star_score"]<bs_max_score)]
# else:
#    filtered_interactions=test_interactions[(test_interactions["star_score"]>bs_min_score) & (test_interactions["star_score"]<bs_max_score)].drop_duplicates(subset="star_score")


