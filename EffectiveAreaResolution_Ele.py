# Effective area from lateral resolution for Electron (W/m2K to nW/K)
# Last modified: July 15 2021
# input: gap [nm], HTC [W/m2K], current density [A/m2]
# Coded by Takuro TOKUNAGA

import numpy as np
import time
import pandas as pd
import sys

sys.path.append('../derjaguin/')
from lip_Ele import liner_interpolate # import function

start = time.time()

# file open
f0 = open('TotalConductance_resolution_interpolate.txt', 'w') # write mode
f1 = open('gapmin.txt', 'r') # read mode
f2 = open('gapmax.txt', 'r') # read mode
f3 = open('radius.txt', 'r') # read mode
f4 = open('TotalConductance_resolution_analytical.txt', 'w') # write mode

# read inputs
# gapmin [m]
for line1 in f1:
    #print(str(line1))
    gapmin = float(line1) # minimum gap distance (input)

# gapmax [m]
for line2 in f2:
    #print(str(line2))
    gapmax = float(line2) # minimum gap distance (input)

# R: Shere radius [m]
for line3 in f3:
    #print(str(line3))
    cr = float(line3) # R (capital r), radius of AFM tip, l:large (input)

# unit conversion
ucnano = 1.0*np.power(10.,-9)
ucangs = 1.0*np.power(10.,-10)

# gap (interatomic distance)
gap = gapmin # initialization, m


# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

def eafr(arg_r, arg_d): # eafr: Effective area from resolution
    resolution_x = np.sqrt(2.0*ucangs*(arg_r+arg_d)) # [m]
    effective_area = np.pi*np.power(resolution_x,2.0)
    cross_section = np.pi*np.power(arg_r,2.0)

    print(arg_r)
    print(effective_area)
    print(cross_section)

    return effective_area

# main
begin()

# parameters
radius_au = 1.740*ucangs # m, atomic radius of Au

#
number = 17 # number of data points x & fx, from 0 to 11
sx=np.zeros((number), dtype='float64')
fx=np.zeros((number), dtype='float64')
fx2=np.zeros((number), dtype='float64')

data = pd.read_csv("../derjaguin/input_h_ele_au.txt", sep="\t", header=None) # in case of electron by WKB
data.columns = ["gap (nm)", "HTC (W/m2K)", "Current Density (A/m2)"]

# input data into tables
for i in range(0, number):
    sx[i] = data.iat[i,0]  # x line, gap [nm]
    fx[i] = data.iat[i,1]  # y1 line, HTC [W/m2K]
    fx2[i] = data.iat[i,2] # y2 line, Current Density [A/m2]

# analytical loop
# analytical solution, not interporation
for i in range(0, number):
    conductance_analytical = fx[i]*eafr(radius_au, sx[i]*ucnano) # [W/K], (W/m2K)*m2
    current_analytical = fx2[i]*eafr(radius_au, sx[i]*ucnano) # [A], (A/m2)*m2

    # Output result
    f4.write(str(sx[i])) # [nm]
    f4.write(str(' '))
    f4.write(str(conductance_analytical/ucnano)) # total conductance [nW/K]
    f4.write(str(' '))
    f4.write(str(current_analytical/ucnano)) # total current [nA]
    f4.write('\n')

# loop for gap
# conductance, current
while gap < gapmax: # [m]
    # initialization
    ind_h = 0
    ind_g = 0
    ind_c = 0  # total conductance
    ind_cd = 0 # total current

    ## thermal conductance [nW/K] ##
    # individual heat transfer coefficient
    ind_h = liner_interpolate(gap/ucnano)[0] # ([nm]), [W/m2K]
    # individual conductance
    ind_g = ind_h*eafr(cr, gap) # [W/K], (W/m2K)*m2

    ## Current [nA] ##
    # individual current density
    ind_cd = liner_interpolate(gap/ucnano)[1] # ([nm]), [A/m2]
    # individual current
    ind_c = ind_cd*eafr(cr, gap) # [A], (A/m2)*m2

    # Output result
    f0.write(str(gap/ucnano)) # [nm]
    f0.write(str(' '))
    f0.write(str(ind_g/ucnano)) # total conductance [nW/K]
    f0.write(str(' '))
    f0.write(str(ind_c/ucnano)) # total current [nA]
    f0.write('\n')

    # current gap display
    print("gap:{:.2f}".format(gap/ucnano) + "[nm]")

    # dgap update (shorter version)
    if gap < 0.098*ucnano: # ~ 0.1 nm
        dgap = 0.01*ucnano # m
    #elif gap>=0.099*ucnano and gap < 0.98*ucnano: # 0.1 nm ~ 1.0 nm
    elif gap>=0.099*ucnano and gap < 1.98*ucnano: # 0.1 nm ~ 1.0 nm
        dgap = 0.1*ucnano # m
    #elif gap>0.99*ucnano and gap < 10.0*ucnano: # 1.0 nm ~ 10.0 nm
    elif gap>1.99*ucnano and gap < 10.0*ucnano: # 1.0 nm ~ 10.0 nm
        dgap = 1.0*ucnano # m
    elif gap>=10*ucnano: # 10.0 nm ~
        dgap = 10*ucnano # m

    # gap update
    gap = gap + dgap

# file close
f0.close()
f1.close()
f2.close()
f3.close()
f4.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
