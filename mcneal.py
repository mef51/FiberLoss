#!/usr/bin/python

# This code reproduces the results from Mcneal's (1976) paper, with ionic current equations
# from Hodgkin and Huxley, based on Boucher's (2011) paper

from __future__ import division
from math import exp, pi
import numpy as np
import pylab

voltShift = -65.5 # mV
Constants = {
    "gBarNa": 120.0, # mS/cm^2
    "gBarK" : 36.0, # mS/cm^2
    "gBarL" : 0.25, # mS/cm^2
    "E_Na"  : 50 - voltShift, # mV
    "E_K"   : -77 - voltShift, # mV
    "E_L"   : -54.4 - voltShift  # mV
}

def alphaN(v):
    if v == 35: return 0.2
    a = 1 - exp((35 - v)/10.0)
    a = a**-1
    a *= 0.02*(v-35)
    return a
alphaN = np.vectorize(alphaN)

def betaN(v):
    if v == 10: return 0.5
    b = (1 - exp((v - 10)/10.0)) ** -1
    b *= 0.05*(10 - v)
    return b
betaN = np.vectorize(betaN)

###### Sodium (Na) Channel (activating)
def alphaM(v):
    if v == 22: return 1.08
    a = (1 - exp((22 - v)/3.0)) ** -1
    a *= 0.36*(v - 22)
    return a
alphaM = np.vectorize(alphaM)

def betaM(v):
    if v == 13: return 8.0
    b = (1 - exp((v - 13)/20.0)) ** -1
    b *= 0.4*(13 - v)
    return b
betaM = np.vectorize(betaM)

###### Sodium (Na) Channel (inactivating)
def alphaH(v):
    if v == -10: return 0.6
    a = (1 - exp((v + 10)/6.0)) ** -1
    a *= 0.1 * (-10 - v)
    return a
alphaH = np.vectorize(alphaH)

def betaH(v):
    b = (1 + exp((45 - v)/10.0)) ** -1
    b *= 4.5
    return b
betaH = np.vectorize(betaH)

def sodiumCurrent(m, h, v):
    # EFRT = (E*F) / (R*T)
    # return pBarNa * h * m**2 * EFRT * F * (NaO - NaI*exp(EFRT)) / (1 - exp(EFRT))
    return m**3 * h * Constants["gBarNa"] * (v - Constants["E_Na"])

def potassiumCurrent(n, v):
    return n**4 * Constants["gBarK"] * (v - Constants["E_K"])

def leakageCurrent(v):
    return Constants["gBarL"] * (v - Constants["E_L"])


nInf  = lambda v: alphaN(v)/(alphaN(v) + betaN(v))
mInf  = lambda v: alphaM(v)/(alphaM(v) + betaM(v))
hInf  = lambda v: alphaH(v)/(alphaH(v) + betaH(v))

def plotAlphaBetaFunctions():
    v = np.arange(-75, 125) # millivolts
    pylab.figure()
    pylab.xlim([-75, 125])
    pylab.plot(v, alphaM(v), v, alphaH(v), v, alphaN(v), v, betaM(v), v, betaH(v), v, betaN(v))
    pylab.legend(('alphaM', 'alphaH', 'alphaN', 'betaM', 'betaH', 'betaN'))
    pylab.title('Alpha and Beta Functions')
    pylab.ylabel(u'Rate Constant (ms^-1)')
    pylab.xlabel('Voltage (mV)')
    pylab.show()

def last(list):
    return list[len(list) - 1]

restingVoltage = 0.0 # mV
dt = 0.025 # ms
T  = 1.0 # ms
cm = 0.0002 # mF/cm^2
D = 0.002 # cm (20microns)
d = 0.7 * D # cm
l = 0.00025 # cm (2.5 microns)

m = [mInf(restingVoltage)]
h = [hInf(restingVoltage)]
n = [nInf(restingVoltage)]

Vm = [restingVoltage]

###############
def printStatus(t, v):
    print "t: " + str(t) + " V: " + str(v)


printStatus(0.0, last(Vm))
for i in range(0, int(T/dt)):
    iNa = sodiumCurrent(last(m), last(h), last(Vm))
    iK  = potassiumCurrent(last(n), last(Vm))
    iL  = leakageCurrent(last(Vm))
    ionicCurrent = iNa + iK + iL

    mNew = (alphaM(last(Vm))*(1 - last(m)) - betaM(last(Vm)) * last(m)) * dt + last(m)
    hNew = (alphaH(last(Vm))*(1 - last(h)) - betaH(last(Vm)) * last(h)) * dt + last(h)
    nNew = (alphaN(last(Vm))*(1 - last(n)) - betaN(last(Vm)) * last(n)) * dt + last(n)
    m.append(mNew)
    h.append(hNew)
    n.append(nNew)

    newV = dt / (cm*pi*d*l) * (0.3 - pi*d*l *(ionicCurrent)) + last(Vm)
    printStatus((i+1)*dt, newV)
    Vm.append(newV)
