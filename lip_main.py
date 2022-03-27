# Calling Liner InterPolation (LIP)
# Last modified: February 01 2019
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys

sys.path.append('../derjaguin/')
from lip import liner_interpolate # import function

start = time.time()

# file open
f0 = open('interpolated_HTC.txt', 'w') # write mode
f1 = open('gapmin.txt', 'r') # read mode
f2 = open('gapmax.txt', 'r') # read mode

# read inputs
# gapmin [m]
for line1 in f1:
    #print(str(line1))
    gapmin = float(line1) # minimum gap distance (input)

# gapmax [m]
for line2 in f2:
    #print(str(line2))
    gapmax = float(line2) # minimum gap distance (input)

# unit conversion
ucnano = 1.0*np.power(10.,-9)

# gap (interatomic distance)
gap = gapmin # initialization, m

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# loop for gap
while gap < gapmax:

    # calling liner interpolation
    ind_h = liner_interpolate(gap/ucnano) # ([nm]), [W/m2K]

    # Output result
    f0.write(str(gap/ucnano)) # [nm]
    f0.write(str(' '))
    f0.write(str(ind_h)) # HTC [W/m2K]
    f0.write('\n')

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

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
