# Derjaguin approximation
# Phys. Rev. B. 94 045406 (2016)
# Last modified: June 11 2018
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
from scipy.integrate import trapz, simps, quad, quadrature, romberg

start = time.time()

sys.path.append('../distance/')
from fe import fe_d   # import function
from agf import agf_d   # import function

# unit conversion
ucnano = 1.0*np.power(10.,-9)

# file open
f = open('TotalHTC.txt', 'w') # write mode
f1 = open('gapmin.txt', 'r') # read mode
f2 = open('gapmax.txt', 'r') # read mode
f3 = open('area.txt', 'r') # read mode

# read inputs
# gapmin
for line1 in f1:
    #print(str(line1))
    gapmin = float(line1) # minimum gap distance (input)

# initialization of gap
gap = gapmin

# gapmax
for line2 in f2:
    #print(str(line2))
    gapmax = float(line2) # minimum gap distance (input)

# area, area of AFM tip
for line3 in f3:
    #print(str(line3))
    area = float(line3) # R, radius of AFM tip, l:large (input)

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# loop for gapmin
while gapmin < gapmax:

    # calculation of h[d,T] for d tilde, call as external functions
    #htc = fe_d(gap) # [w/m2K]
    htc = agf_d(gap) # [w/m2K]

    gc = area*htc # conductance # [w/K]

    # Output result
    f.write(str(gap/ucnano)) # [nm]
    f.write(str(' '))
    f.write(str(gc)) # total conductance [W/K]
    f.write('\n')

    # update of dgap (gap distance)
    if gap < 0.1*ucnano:
        dgap = 0.01*ucnano
    elif gap>=0.1*ucnano and gap < 1*ucnano:
        dgap = 0.1*ucnano
    elif gap>=1*ucnano and gap < 10*ucnano:
        dgap = 1*ucnano
    elif gap>=10*ucnano:
        dgap = 10*ucnano

    # current gap
    print("gap:{:.2f}".format(gap/ucnano) + "[nm]")

    # update of gap
    gap = gap + dgap

# file close
f.close()
f1.close()
f2.close()
f3.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
