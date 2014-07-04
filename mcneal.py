#!/usr/bin/python

# This code reproduces the results from Mcneal's (1976) paper, with ionic current equations
# from Hodgkin and Huxley, based on Boucher's (2011) paper

from __future__ import division
from scipy.integrate import odeint
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
    return 0.02*(v-35) * ((1 - exp((35 - v)/10.0)) ** -1)
alphaN = np.vectorize(alphaN)

def betaN(v):
    if v == 10: return 0.5
    return 0.05*(10 - v) * ((1 - exp((v - 10)/10.0)) ** -1)
betaN = np.vectorize(betaN)

###### Sodium (Na) Channel (activating)
def alphaM(v):
    if v == 22: return 1.08
    return 0.36*(v - 22) * ((1 - exp((22 - v)/3.0)) ** -1)
alphaM = np.vectorize(alphaM)

def betaM(v):
    if v == 13: return 8.0
    return 0.4*(13 - v) * ((1 - exp((v - 13)/20.0)) ** -1)
betaM = np.vectorize(betaM)

###### Sodium (Na) Channel (inactivating)
def alphaH(v):
    if v == -10: return 0.6
    return 0.1 * (-10 - v) * ((1 - exp((v + 10)/6.0)) ** -1)
alphaH = np.vectorize(alphaH)

def betaH(v):
    return 4.5 * ((1 + exp((45 - v)/10.0)) ** -1)
betaH = np.vectorize(betaH)

###### Non-specific Delayed Current Density (sodium!)
def alphaP(v):
    if v == 40: return 0.06
    return 0.006 * (v-40) * ((1 - exp((40 - v)/10)) ** -1)
alphaP = np.vectorize(alphaP)

def betaP(v):
    if v == -25: return 1.8
    return 0.09 * (-25 - v) * ((1 - exp((v + 25)/20)) ** -1)
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

def iNa(m, h, v):
    E = (v + Constants["Vr"])/1000; F = Constants["F"]; R = Constants["R"]; T = Constants["T"]
    NaO = Constants["NaO"]  # (mol/cm^3)
    NaI = Constants["NaI"]  # (mol/cm^3)
    pBarNa = 8e-3 # cm/s
    EFRT = (E*F) / (R*T)
    return pBarNa * h * m**2 * EFRT * F * (NaO - NaI*exp(EFRT)) / (1 - exp(EFRT))

def iK(n, v):
    E = (v + Constants["Vr"])/1000; F = Constants["F"]; R = Constants["R"]; T = Constants["T"]
    Ko = 2.5 / 1e6 # (mol/cm^3)
    Ki = 120 / 1e6 # (mol/cm^3)
    pBarK = 1.2e-3 # cm/s
    EFRT = (E*F) / (R*T)
    return pBarK * n**2 * EFRT * F * (Ko - Ki*exp(EFRT)) / (1 - exp(EFRT))

def iL(v):
    return Constants["gBarL"] * (v - Constants["E_L"])

def iP(p, v):
    E = (v + Constants["Vr"])/1000; F = Constants["F"]; R = Constants["R"]; T = Constants["T"]
    NaO = Constants["NaO"]  # (mol/cm^3)
    NaI = Constants["NaI"]  # (mol/cm^3)
    pBarP = 0.54e-3 # cm/s
    EFRT = (E*F) / (R*T)
    return pBarP * p**2 * EFRT * F * (NaO - NaI*exp(EFRT)) / (1 - exp(EFRT))

def iIonic(m, h, n, p, v):
    return iNa(m, h, v) + iK(n, v) + iL(v) + iP(p, v)

nInf  = lambda v: alphaN(v)/(alphaN(v) + betaN(v))
mInf  = lambda v: alphaM(v)/(alphaM(v) + betaM(v))
hInf  = lambda v: alphaH(v)/(alphaH(v) + betaH(v))
pInf  = lambda v: alphaP(v)/(alphaP(v) + betaP(v))

restingVoltage = -0.1888767 # mV
dt = 0.00025 # ms
T  = 1 # ms
cm = 0.002 # mF/cm^2
D = 0.002 # cm (20microns)
d = 0.7 * D # cm
l = 0.00025 # cm (2.5 microns)
r = 0.1  # cm (1mm)
I = 0.3 # mA
RhoE = 300e3 # mohm*cm
RhoI = 110e3 # mohm*cm
L = 0.2 # cm
Ga = 1 / (RhoI * L)
Ve = RhoE * I / (4 * pi * r)

m = mInf(restingVoltage)
h = hInf(restingVoltage)
n = nInf(restingVoltage)
p = pInf(restingVoltage)
Vm = restingVoltage

initialValues = [Vm, m, h, n, p]
timeLine = [i*dt for i in range(0, int(T/dt))]
abserr = 1.0e-7
relerr = 1.0e-7

def f(values, t):
    v, m, h, n, p = values
    funcs = [
        (1 / (cm)) * (-2*Ga*(Ve + v) - iIonic(m, h, n, p, v)),
        alphaM(v)*(1 - m) - betaM(v) * m,
        alphaH(v)*(1 - h) - betaH(v) * h,
        alphaN(v)*(1 - n) - betaN(v) * n,
        alphaP(v)*(1 - p) - betaP(v) * p
    ]
    return funcs

solution = odeint(f, initialValues, timeLine, atol = abserr, rtol=relerr)

pylab.plot(timeLine, solution[:,0])
pylab.ylabel('V (mV)')
pylab.xlabel('Time, t (ms)')
pylab.show()
