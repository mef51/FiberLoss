#!/usr/bin/python

# This code reproduces the results from Mcneal's (1976) paper, with ionic current equations
# from Hodgkin and Huxley, based on Boucher's (2011) paper

from __future__ import division
from math import exp, pi
import numpy as np
import pylab

Constants = {
    "Vr"    : -70,     # mV
    "F"     : 964853.0,# C/mole
    "R"     : 8.3144,  # J/K/mole
    "T"     : 295.18,  # K
    "gBarL" : 30.3, # mS/cm^2
    "E_L"   : 0.026, # mV
    "NaO"   : 114.5 / 1e6, # (mol/cm^3)
    "NaI"   : 13.7 / 1e6  # (mol/cm^3)
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

###### Non-specific Delayed Current Density (sodium!)
def alphaP(v):
    if v == 40: return 0.06
    a = (1 - exp((40 - v)/10)) ** -1
    a *= 0.006 * (v-40)
    return a
alphaP = np.vectorize(alphaP)

def betaP(v):
    if v == -25: return 1.8
    b = (1 - exp((v + 25)/20)) ** -1
    b *= 0.09 * (-25 - v)
    return b
betaP = np.vectorize(betaP)

def plotAlphaBetaFunctions():
    v = np.arange(-75, 125) # millivolts
    pylab.figure()
    pylab.xlim([-75, 125])
    pylab.plot(v, alphaM(v), v, alphaH(v), v, alphaN(v), v, alphaP(v), '--', v, betaM(v), v, betaH(v), v, betaN(v), v, betaP(v), '--')
    pylab.legend(('alphaM', 'alphaH', 'alphaN', 'alphaP', 'betaM', 'betaH', 'betaN', 'betaP'))
    pylab.title('Forward and Backward Rates')
    pylab.ylabel(u'Rate Constant (ms^-1)')
    pylab.xlabel('Voltage (mV)')
    pylab.show()

def sodiumCurrent(m, h, v):
    E = v + Constants["Vr"]; F = Constants["F"]; R = Constants["R"]; T = Constants["T"]
    NaO = Constants["NaO"]  # (mol/cm^3)
    NaI = Constants["NaI"]  # (mol/cm^3)
    pBarNa = 8e-3 # cm/s
    EFRT = (E*F) / (R*T)
    return pBarNa * h * m**2 * EFRT * F * (NaO - NaI*exp(EFRT)) / (1 - exp(EFRT))

def potassiumCurrent(n, v):
    E = v + Constants["Vr"]; F = Constants["F"]; R = Constants["R"]; T = Constants["T"]
    Ko = 2.5 / 1e6 # (mol/cm^3)
    Ki = 120 / 1e6 # (mol/cm^3)
    pBarK = 1.2e-3 # cm/s
    EFRT = (E*F) / (R*T)
    return pBarK * n**2 * EFRT * F * (Ko - Ki*exp(EFRT)) / (1 - exp(EFRT))

def leakageCurrent(v):
    return Constants["gBarL"] * (v - Constants["E_L"])

def delayedCurrent(p, v):
    E = v + Constants["Vr"]; F = Constants["F"]; R = Constants["R"]; T = Constants["T"]
    NaO = Constants["NaO"]  # (mol/cm^3)
    NaI = Constants["NaI"]  # (mol/cm^3)
    pBarP = 0.54e-3 # cm/s
    EFRT = (E*F) / (R*T)
    return pBarP * p**2 * EFRT * F * (NaO - NaI*exp(EFRT)) / (1 - exp(EFRT))

nInf  = lambda v: alphaN(v)/(alphaN(v) + betaN(v))
mInf  = lambda v: alphaM(v)/(alphaM(v) + betaM(v))
hInf  = lambda v: alphaH(v)/(alphaH(v) + betaH(v))
pInf  = lambda v: alphaP(v)/(alphaP(v) + betaP(v))

def last(list):
    return list[len(list) - 1]

restingVoltage = 0.0 # mV
dt = 0.025 # ms
T  = 1.0 # ms
cm = 0.0002 # mF/cm^2
D = 0.002 # cm (20microns)
d = 0.7 * D # cm
l = 0.00025 # cm (2.5 microns)
r = 0.1  # cm (1mm)
I = 0.3  # mA
RhoE = 300 # ohm*cm
RhoI = 110 # ohm*cm
L = 0.2 # cm

m = [mInf(restingVoltage)]
h = [hInf(restingVoltage)]
n = [nInf(restingVoltage)]
p = [pInf(restingVoltage)]

Vm = [restingVoltage]

###############
def printStatus(t, v):
    print "t: " + str(t) + " V: " + str(v)

def printVar(v, name):
    print name + ": " + str(v)

printStatus(0.0, last(Vm))
Ga = 1 / (RhoI * L)
for i in range(0, int(T/dt)):
    # printVar(m, 'm')
    # printVar(h, 'h')
    # printVar(n, 'n')
    # printVar(p, 'p')
    # printVar(Vm, 'Vm')

    iNa = sodiumCurrent(last(m), last(h), last(Vm))
    iK  = potassiumCurrent(last(n), last(Vm))
    iL  = leakageCurrent(last(Vm))
    iP  = delayedCurrent(last(p), last(Vm))
    ionicCurrent = iNa + iK + iL + iP

    mNew = (alphaM(last(Vm))*(1 - last(m)) - betaM(last(Vm)) * last(m)) * dt + last(m)
    hNew = (alphaH(last(Vm))*(1 - last(h)) - betaH(last(Vm)) * last(h)) * dt + last(h)
    nNew = (alphaN(last(Vm))*(1 - last(n)) - betaN(last(Vm)) * last(n)) * dt + last(n)
    pNew = (alphaP(last(Vm))*(1 - last(p)) - betaP(last(Vm)) * last(p)) * dt + last(p)
    m.append(mNew)
    h.append(hNew)
    n.append(nNew)
    p.append(pNew)

    Ve = RhoE * I / (4 * pi * r)

    print iNa
    print iK
    print iL
    print iP
    newV = dt / (cm) * (-2*Ga*(Ve + last(Vm)) - (ionicCurrent)) + last(Vm)
    printStatus((i+1)*dt, newV)
    Vm.append(newV)
