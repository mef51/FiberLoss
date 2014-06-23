#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import pylab

##########################
# parameters

# distribution of the diameters of the fibers in the nerve
# the units are micrometers (µm) to relative frequency (%)
nerve = {
    1  : 1/17.0,
    2  : 1/17.0,
    3  : 1/17.0,
    4  : 1/17.0,
    5  : 1/17.0,
    6  : 1/17.0,
    7  : 1/17.0,
    8  : 1/17.0,
    9  : 1/17.0,
    10 : 1/17.0,
    11 : 1/17.0,
    12 : 1/17.0,
    13 : 1/17.0,
    14 : 1/17.0,
    15 : 1/17.0,
    16 : 1/17.0,
    17 : 1/17.0
}

totalNumFibers = 187 # some number divisible by len(nerve)

def r(n):
    return 1.0

rho_e = 2
fibreDiameter = 10 # µm

# stuff
axonalDiameter = 0.76 * fibreDiameter + 1.37e-7

# `n`
def numFibers(diameter):
    return nerve[diameter] * totalNumFibers

# `L`
def internodalLength(D):
    if D <= 1.2e-5:
        return 102*D + 7.15e-5
    else:
        return 4.68e-4 * np.log(D / 9.74e-7)

# `Ge`
def externalConductance(n):
    return 1 / (4 * np.pi * r(n) * rho_e)

## Current stuff
def e(D):
    return S(D) / (np.pi * (D/2)**2)

def i_D(x):
    # sum of currents at position x
    return -1

def InD(x):
    numFibers(fibreDiameter) * e(fibreDiameter) * i_D(x)


#############
# Integrate
dt = 0.05 # ms
simulationTime = 10 # ms
timeLine = np.arange(0, simulationTime + dt, dt)
for i, t in enumerate(timeLine):
    print t
