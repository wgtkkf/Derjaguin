# Derjaguin approximation for bumpy surface
# Last modified: July 03 2018
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import pandas as pd
import sys
from scipy.integrate import trapz, simps, quad, quadrature, romberg

sys.path.append('../regression/')
from fitting import fitting # read fitting function

start = time.time()

sys.path.append('../distance/')
from fe import fe_d   # import function
from agf import agf_d   # import function

# unit conversion
ucnano = 1.0*np.power(10.,-9)

# parameters
numx = 17#32
numy = 17#32
numxy = numx*numy

#
sx=np.zeros((numxy), dtype='float64')

# file open
f = open('TotalConductance.txt', 'w') # write mode
f1 = open('gapmin.txt', 'r') # read mode
f2 = open('gapmax.txt', 'r') # read mode
f3 = open('area.txt', 'r') # read mode
f4 = open('datasheet.txt', 'w')

# read inputs
# gapmin
for line1 in f1:
    #print(str(line1))
    gapmin = float(line1) # minimum gap distance (input)

# gapmax
for line2 in f2:
    #print(str(line2))
    gapmax = float(line2) # minimum gap distance (input)

# area, area of AFM tip
for line3 in f3:
    #print(str(line3))
    area = float(line3) # R, radius of AFM tip, l:large (input)

# initialization of gap related parameters
cut_gap = 0.1*np.power(10.,-9) # [nm]
htc = 0

gap = gapmin
gap_ch = gapmin # characteristic gap
dgap_ch = gapmin
gap_limit = 0.5*ucnano # [m]

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

#bumpy_data = pd.read_csv("../fe/text/flux"+str(j)+".txt", sep="\t", header=None) # change file name here tab:\t
bumpy_data = pd.read_csv("../surface/bumpy/data0.txt", sep=" ", header=None)
bumpy_data.columns = ["x coord", "y coord", "bumpy"]

# input data into tables
for i in range(0, numxy):
    sx[i] = bumpy_data.iat[i,2]*ucnano # bumpy from [nm] to [m]
    #print(sx[i])

# loop for gapmin
while gap_ch <= gapmax:
#while gap_ch < gapmax:

    for i in range(0, numxy): # bumpy loop
        gap = gap_ch-sx[i]

        if gap < cut_gap: # below cut_gap, conductance gives same value
            gap = cut_gap

        # consider area factor below gap_limit
        if gap>=gap_limit:
            temp = agf_d(gap)       # [W/m2K]
            htc = htc + temp        # [W/m2K] for each bumpy
        elif gap<gap_limit:         # below gap limit
            temp = agf_d(gap_limit) # [W/m2K]
            htc = htc + temp        # [W/m2K] for each bumpy

        # file output of conductances for each bumpy & fitted conductance
        f4.write(str(temp)) # conductance [W/m2K]
        f4.write(str(' '))

    # fitting value for each gap_ch
    fitted_conductance = fitting(gap_ch/ucnano) # gap_ch: input as [nm], [nW/K]
    f4.write(str(fitted_conductance*ucnano))    # fitted conductance,    [W/K]
    f4.write('\n')

    # averaged total conductance @ each gap_ch
    gc = area*(htc/numxy)                       # conductance # [W/K]

    # Output result
    f.write(str(gap_ch/ucnano)) # [nm]
    f.write(str(' '))
    f.write(str(gc/ucnano)) # total conductance [nW/K]
    f.write('\n')

    # current gap (characteristic gap)
    print("gap_ch:{:.2f}".format(gap_ch/ucnano) + "[nm]")

    # update of dgap_ch (gap distance)
    if gap_ch < 0.09*ucnano:
        dgap_ch = 0.01*ucnano
    elif gap_ch>=0.1*ucnano and gap_ch < 1*ucnano:
        dgap_ch = 0.1*ucnano
    elif gap_ch>=1*ucnano and gap_ch < 10*ucnano:
        dgap_ch = 1*ucnano
    elif gap_ch>=10*ucnano:
        dgap_ch = 10*ucnano

    # gap update
    gap_ch = gap_ch + dgap_ch

    # reset parameters
    htc = 0
    gc = 0
    temp = 0

# file close
f.close()
f1.close()
f2.close()
f3.close()
f4.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
