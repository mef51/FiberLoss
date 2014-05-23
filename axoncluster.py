#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

class AxonPositionNode:
    """A node of Ranvier on an Axon, as modelled by Hodgkin and Huxley in 1952"""
    def __init__(self, x, y):
        # position
        self.x = x
        self.y = y

        # Potassium (K) Channel
        alphaN = self.alphaN = np.vectorize(lambda v: 0.01*(10 - v) / (np.exp((10-v)/10) - 1) if v != 10 else 0.1)
        betaN = self.betaN  = lambda v: 0.125 * np.exp(-v/80)
        nInf = self.nInf   = lambda v: alphaN(v)/(alphaN(v) + betaN(v))

        # Sodium (Na) Channel (activating)
        alphaM = self.alphaM = np.vectorize(lambda v: 0.1*(25-v) / (np.exp((25-v)/10) - 1) if v!= 25 else 1)
        betaM = self.betaM = lambda v: 4 * np.exp(-v/18)
        mInf = self.mInf = lambda v: alphaM(v)/(alphaM(v) + betaM(v))

        # Sodium (Na) Channel (inactivating)
        alphaH = self.alphaH = lambda v: 0.07 * np.exp(-v/20)
        betaH = self.betaH  = lambda v: 1/(np.exp((30-v)/10) + 1)
        hInf = self.hInf   = lambda v: alphaH(v)/(alphaH(v) + betaH(v))

        # Hodgkin-Huxley Parametahs (from the papah!)
        params = self.params = {
            "restingVoltage"     : 0,      # V_rest (mv)
            "Cm"                 : 1,      # uF/cm2
            "gBarNa"             : 120,    # mS/cm2
            "gBarK"              : 36,     # mS/cm2
            "gBarL"              : 0.3,    # mS/cm2
            "sodiumPotential"    : 115,    # mV
            "potassiumPotential" : -12,    # mv
            "leakagePotential"   : 10.613 # mV
        }

        self.Vm    = [params["restingVoltage"]] # The axon node's membrane potential
        self.m     = mInf(params["restingVoltage"])
        self.h     = hInf(params["restingVoltage"])
        self.n     = nInf(params["restingVoltage"])

    # integrate response to stimulus current `stimulus`
    def step(self, stimulus, dt):
        I = stimulus # I[i-1]
        lastVm = len(self.Vm) - 1

        sodiumConductance    = self.params["gBarNa"] * (self.m**3) * self.h
        potassiumConductance = self.params["gBarK"]  * (self.n**4)
        leakageConductance   = self.params["gBarL"]

        # integrate the equations on m, h, and n
        self.m += (self.alphaM(self.Vm[lastVm]) * (1 - self.m) - self.betaM(self.Vm[lastVm])*self.m) * dt
        self.h += (self.alphaH(self.Vm[lastVm]) * (1 - self.h) - self.betaH(self.Vm[lastVm])*self.h) * dt
        self.n += (self.alphaN(self.Vm[lastVm]) * (1 - self.n) - self.betaN(self.Vm[lastVm])*self.n) * dt

        # now integrate the changes in V
        sodiumCurrent = sodiumConductance * (self.Vm[lastVm] - self.params["sodiumPotential"])
        potassiumCurrent = potassiumConductance * (self.Vm[lastVm] - self.params["potassiumPotential"])
        leakageCurrent = leakageConductance * (self.Vm[lastVm] - self.params["leakagePotential"])
        self.Vm.append(self.Vm[lastVm] + (I - sodiumCurrent - potassiumCurrent - leakageCurrent) * dt / self.params["Cm"])


# prototypin
class Simulation:
    def __init__(self, T, dt):
        self.T = T
        selft.dt = dt
        self.timeLine = np.arange(0, T+dt, dt)

T    = 55    # ms
dt   = 0.025 # ms
timeLine = np.arange(0, T+dt, dt)

# Current Stimulus
def getCurrent(t, distance, stimulusCurrent, tPulseStart=5, pulseWidth=25):
    if tPulseStart <= t <= (tPulseStart + pulseWidth):
        return stimulusCurrent / (2*np.pi * distance**2) # uA/cm2. the current that reaches the axon.
    else:
        return 0


numAxons = 10
stimulus = 10 # uA
# distance = 1 # cm

# the nerve is a bundle of nerve fibers, or, a cluster of axons
nerveDistance = 1 # cm. The distance the center of the nerve is from the stimulus
nerveRadius = 2 # cm

# axons = []
# for i in range(0, numAxons):
#     axons.append(AxonNode(distance))
