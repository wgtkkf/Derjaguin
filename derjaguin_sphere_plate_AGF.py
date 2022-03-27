# Derjaguin approximation for AGF (W/m2K to nW/K)
# Nat. Photonics. 3 514 (2009)
# Last modified: April 01 2019
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys

sys.path.append('../derjaguin/')
from lip_AGF import liner_interpolate # import function

start = time.time()

# file open
f0 = open('TotalConductance.txt', 'w') # write mode
f1 = open('gapmin.txt', 'r') # read mode
f2 = open('gapmax.txt', 'r') # read mode
f3 = open('radius.txt', 'r') # read mode

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

# gap (interatomic distance)
gap = gapmin # initialization, m

# parameters for radius
nrmax = 1000 # division nrmax toward r direction, integer type
srmin = 0 # initialization of rmin, s:small
sr = srmin # initialization of r, s:small
srmax = cr # initialization of rmax, s:small
dsr = (srmax-srmin)/nrmax # initialization of dr, s:small

# conductance
total_g = 0 # total conductance

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

def dtilde(r, d):
    y = d + cr - np.sqrt(np.power(cr,2)-np.power(r,2))

    return y

# main
begin()

# loop for gap
while gap < gapmax:
    # initialization of the parameters
    dtldr = 0
    ind_h = 0
    total_g = 0

    # loop for r
    while sr < srmax:
        # d tilder (dtldr)
        dtldr = dtilde(sr, gap) # [m]

        # individual heat transfer coefficient
        ind_h = liner_interpolate(dtldr/ucnano) # ([nm]), [W/m2K]

        # individual conductance
        ind_g = ind_h*(2*np.pi*sr)*dsr # [W/K]

        # conductance (summation)
        total_g = total_g + ind_g # [W/K]

        # sr update
        sr = sr + dsr # [m]

    # Output result
    f0.write(str(gap/ucnano)) # [nm]
    f0.write(str(' '))
    f0.write(str(total_g/ucnano)) # total conductance [nW/K]
    f0.write('\n')

    # reset the parameters
    sr = srmin
    total_g = 0

    # current gap display
    print("gap:{:.2f}".format(gap/ucnano) + "[nm]")

    # dgap update (shorter version)
    if gap < 0.098*ucnano: # ~ 0.1 nm
        dgap = 0.01*ucnano # m
    elif gap>=0.099*ucnano and gap < 0.98*ucnano: # 0.1 nm ~ 1.0 nm
        dgap = 0.1*ucnano # m
    elif gap>0.99*ucnano and gap < 10.0*ucnano: # 1.0 nm ~ 10.0 nm
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

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
