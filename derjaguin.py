# Derjaguin approximation
# Nat. Photonics. 3 514 (2009)
# Last modified: June 06 2018
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

# division nrmax toward r direction
nrmax = 100 # use int number

# file open
f = open('TotalHTC.txt', 'w') # write mode
f1 = open('gapmin.txt', 'r') # read mode
f2 = open('gapmax.txt', 'r') # read mode
f3 = open('lr.txt', 'r') # read mode

# read inputs
# gapmin
for line1 in f1:
    #print(str(line1))
    gapmin = float(line1) # minimum gap distance (input)

# gapmax
for line2 in f2:
    #print(str(line2))
    gapmax = float(line2) # minimum gap distance (input)

# R: Shere radius
for line3 in f3:
    #print(str(line3))
    lr = float(line3) # R, radius of AFM tip, l:large (input)

# discretization points
table_sr = np.zeros(nrmax,dtype='float64') # radius
table_h = np.zeros(nrmax,dtype='float64') # HTC

#  matrix and vectors for 3 order spline
ma=np.zeros((nrmax,nrmax), dtype='float64')
vx=np.zeros((nrmax), dtype='float64')
vb=np.zeros((nrmax), dtype='float64')

# parameters
srmin = 0 # initialization of rmin, s:small
sr = srmin # initialization of sr, s:small
srmax = lr
dsr = (srmax-srmin)/nrmax
total_h = 0 # total conductance

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

def dtilde(r, dmin):
    y = dmin + lr - np.sqrt(np.power(lr,2)-np.power(r,2))

    return y

# main
begin()

# loop for gapmin
while gapmin < gapmax:

    # loop for r
    for j in range(0,nrmax):
        # calculation of d tilde
        dtld = dtilde(sr, gapmin) # d tilde

        if dtld<0:
            print(str('negative d'))
            sys.exit()

        # calculation of h[d,T] for d tilde, call as external functions
        htc = fe_d(dtld)
        #htc = agf_d(dtld)

        # calculation of 'h[d,T]*(2*pi*r)' equation (2)
        table_sr[j] = sr
        table_h[j] = htc*(2*np.pi*sr) # table h

        # distance update
        sr = sr + dsr

    ## integral calculation by 3 order spline
    # initialization of matrix & vectors
    # matrix A
    for i in range(0, nrmax): # from 0 to 13, total 14
        # sub-diagonal components
        if i>0 and i<nrmax-1: # 1 to 12, nrmax=14
            ma[i][i-1] = table_sr[i] - table_sr[i-1] # left
            ma[i][i+1] = table_sr[i+1] - table_sr[i] # right

        #diagonal components
        if i==0:
            ma[i][i] = 1
        elif i>0 and i<nrmax-1: # 1 to 12
            ma[i][i] = 2*(ma[i][i-1]+ma[i][i+1])
        elif i==nrmax-1: # i=13
            ma[i][i] = 1

    # vector b
    for i in range(0, nrmax): # from 0 to 13, total 14
        if i>0 and i<nrmax-1: # from 1 to 12
            h1 = table_sr[i+1]-table_sr[i]
            h0 = table_sr[i]-table_sr[i-1]
            a2 = table_h[i+1]
            a1 = table_h[i]
            a0 = table_h[i-1]

            vb[i] = (3/h1)*(a2-a1)-(3/h0)*(a1-a0)

        # reset parameters
        h1 = 0
        h0 = 0
        a2 = 0
        a1 = 0
        a0 = 0

    # calculation of vector x, LU method
    vx = np.linalg.solve(ma, vb)

    # integral calculation
    for i in range(0, nrmax-1):

        # reset coefficients
        ipfx = 0
        sa = 0
        sh = 0
        sb = 0
        sd = 0

        # interporation function
        def ipfunction(x):
            tempx1 = table_sr[i]
            tempx2 = table_sr[i+1]

            # coefficients
            sa = table_h[i] # small a
            sh = table_sr[i+1]-table_sr[i] # small h
            sb = (1/sh)*(table_h[i+1]-table_h[i]) - (sh/3)*(2*vx[i]+vx[i+1]) # small b
            sc = vx[i] # small c
            sd = (vx[i+1]-vx[i])/(3*sh) # small d

            # define of interporation function
            ipfx = sa + sb*(x-table_sr[i]) + sc*np.power((x-table_sr[i]),2)\
            + sd*np.power((x-table_sr[i]),3)

            return ipfx

        integral = quad(ipfunction, table_sr[i], table_sr[i+1])[0]
        total_h = total_h + integral

    # Output result
    f.write(str(gapmin/ucnano))
    f.write(str(' '))
    f.write(str(total_h)) # total conductance
    f.write('\n')

    # reset all the parameters
    dtld = 0
    sr = 0
    total_h = 0

    # current gap display
    print("gap:{:.2f}".format(gapmin/ucnano) + "[nm]")

    # update of gapmin (gap distance)
    if gapmin < 0.1*ucnano:
        gapmin = 0.1*ucnano
    elif gapmin>=0.1*ucnano and gapmin < 1*ucnano:
        gapmin = 1.0*ucnano
    elif gapmin>=1*ucnano and gapmin < 10*ucnano:
        gapmin = 10*ucnano
    elif gapmin>=10*ucnano:
        gapmin = gapmin + 10*ucnano # +10 [nm]

# file close
f.close()
f1.close()
f2.close()
f3.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
