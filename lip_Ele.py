# Liner InterPolation: LIP for Electron tunneling & Thermionic
# read gap & HTC W/m2K files as input
# Coded by Takuro TOKUNAGA
# Last modified: October 08 2019

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
number = 17 # number of data points x & fx, from 0 to 11
lc_au = 407.82*ucpico # m, lattice constant of Au
cutoff = lc_au/ucnano # [nm]
cutoff = 0.40 # [nm]

#
sx=np.zeros((number), dtype='float64')
fx=np.zeros((number), dtype='float64')
fx2=np.zeros((number), dtype='float64')

# main
#begin()

data = pd.read_csv("../derjaguin/input_h_ele_au.txt", sep="\t", header=None) # in case of electron by WKB
data.columns = ["gap (nm)", "HTC (W/m2K)", "Current Density (A/m2)"]

# input data into tables
for i in range(0, number):
    sx[i] = data.iat[i,0] # x line
    fx[i] = data.iat[i,1] # y1 line
    fx2[i] = data.iat[i,2] # y2 line

# interporation function
def liner_interpolate(arg_x): # arg_x: [nm]
    ipfx = 0
    ipfx2 = 0

    for i in range(0, number-1): # 74 of 0 to 73
        if arg_x < cutoff: # gap < cutoff [nm]
            ipfx = fx[0]
            ipfx2 = fx2[0]

        elif arg_x >= sx[i] and arg_x <=sx[i+1]: # cutoff [nm] < gap < 100 [nm]
            # HTC
            ipfx = fx[i] + (fx[i+1]-fx[i])*(arg_x-sx[i])/(sx[i+1]-sx[i])

            # Current density
            ipfx2 = fx2[i] + (fx2[i+1]-fx2[i])*(arg_x-sx[i])/(sx[i+1]-sx[i])

        elif arg_x > sx[number-1]: # gap > 100 [nm]
            ipfx = 0
            ipfx2 = 0

    return ipfx, ipfx2

#end()

# time display
elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
