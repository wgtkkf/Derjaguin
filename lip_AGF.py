# Liner InterPolation: LIP
# read gap & flux data files as input
# as external function
# Coded by Takuro TOKUNAGA
# Last modified: July 22 2019

import math
import numpy as np
import cmath
import time
import pandas as pd
start = time.time()

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# Unit conversion:
ucev = 1.602176620898*np.power(10.,-19)
ucnano = 1.0*np.power(10.,-9)
ucangs = 1.0*np.power(10.,-10)
ucpico = 1.0*np.power(10.,-12)

# parameters
number = 76#37 # number of data points x & fx, from 0 to 74
lc_au = 407.82*ucpico # m, lattice constant of Au
cutoff = lc_au/ucnano # [nm]

#
sx=np.zeros((number), dtype='float64')
fx=np.zeros((number), dtype='float64')

# main
#begin()

data = pd.read_csv("../derjaguin/input_h_AGF_au.txt", sep="\t", header=None) # in case of AGF
data.columns = ["gap", "heat transfer coefficient"]

# input data into tables
for i in range(0, number):
    sx[i] = data.iat[i,0] # x line
    fx[i] = data.iat[i,1] # fx line

# interporation function
def liner_interpolate(arg_x): # arg_x: [nm]
    ipfx = 0

    for i in range(0, number-1): # 74 of 0 to 73
        if arg_x < cutoff: # gap < cutoff [nm]
            ipfx = fx[0]

        elif arg_x >= sx[i] and arg_x <=sx[i+1]: # cutoff [nm] < gap < 100 [nm]
            ipfx = fx[i] + (fx[i+1]-fx[i])*(arg_x-sx[i])/(sx[i+1]-sx[i])

        elif arg_x > sx[number-1]: # gap > 100 [nm]
            ipfx = 0

    return ipfx

#end()

# time display
elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
