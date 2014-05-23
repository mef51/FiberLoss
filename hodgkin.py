#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pylab
import time
from matplotlib.animation import ArtistAnimation

# sauce: http://www.neurdon.com/2011/01/26/neural-modeling-with-python-part-2/

class AxonNode:
    """A node of Ranvier on an Axon, as modelled by Hodgkin and Huxley in 1952"""
    def __init__(self, distance):

        self.distance = distance # distance from stimulus current

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
    def step(self, stimulus):
        I = stimulus # I[i-1]
        sodiumConductance    = self.params["gBarNa"] * (self.m**3) * self.h
        potassiumConductance = self.params["gBarK"]  * (self.n**4)
        leakageConductance   = self.params["gBarL"]

        # integrate the equations on m, h, and n
        self.m += (self.alphaM(self.Vm[i-1]) * (1 - self.m) - self.betaM(self.Vm[i-1])*self.m) * dt
        self.h += (self.alphaH(self.Vm[i-1]) * (1 - self.h) - self.betaH(self.Vm[i-1])*self.h) * dt
        self.n += (self.alphaN(self.Vm[i-1]) * (1 - self.n) - self.betaN(self.Vm[i-1])*self.n) * dt

        # now integrate the changes in V
        sodiumCurrent = sodiumConductance * (self.Vm[i-1] - self.params["sodiumPotential"])
        potassiumCurrent = potassiumConductance * (self.Vm[i-1] - self.params["potassiumPotential"])
        leakageCurrent = leakageConductance * (self.Vm[i-1] - self.params["leakagePotential"])
        self.Vm.append(self.Vm[i-1] + (I - sodiumCurrent - potassiumCurrent - leakageCurrent) * dt / self.params["Cm"])

# Channel Activity
# v = np.arange(-50, 151) # millivolts
# pylab.figure()
# pylab.plot(v, mInf(v), v, hInf(v), v, nInf(v))
# pylab.legend(('m', 'h', 'n'))
# pylab.title('Steady state values of ion channel gating variables')
# pylab.ylabel('Magnitude')
# pylab.xlabel('Voltage (mV)')
# pylab.savefig("mhn.jpg")

# Setup parameters and state variables
T    = 55    # ms
dt   = 0.025 # ms
timeLine = np.arange(0, T+dt, dt)

# Current Stimulus
def getCurrent(t, distance, stimulusCurrent, tPulseStart=5, pulseWidth=25):
    if tPulseStart <= t <= (tPulseStart + pulseWidth):
        return stimulusCurrent / (2*np.pi * distance**2) # uA/cm2. the current that reaches the axon.
    else:
        return 0

# Main loop
animationFigure = pylab.figure(figsize=(16, 9))
images = []
axons = []
stimulusCurrent = 10. # uA. the current applied at the surface

# create the axons for each plot
for i in range(0, 9):
    distance = (i+1)*0.1 # cm. The distance between the axon and the current
    axons.append(AxonNode(distance))

for i in range(1, len(timeLine)):
    lines = []
    for k in range(1, 10):
        pylab.subplot(3, 3, k)
        axon = axons[k-1]
        # the current at time NOW
        I = getCurrent(i*dt, axon.distance, stimulusCurrent)

        # for the plot's title
        effectiveCurrent = stimulusCurrent / (2*np.pi * axon.distance**2) # uA/cm2. the current that reaches the axon.
        axon.step(I)

        # update status every 25 frames (about a second of video).
        if i % 25 == 0:
            # get a profile of the current up til now
            curr = []
            for j in range(0, len(timeLine[:i+1])):
                curr.append(getCurrent(j*dt, axon.distance, stimulusCurrent))

            # plot a frame of the graph
            voltageLine, currentLine = pylab.plot(timeLine[:i+1], axon.Vm[:i+1], 'b-', timeLine[:i+1], curr[:i+1], 'g-')
            pylab.legend([voltageLine, currentLine], ["Response from cell", "Impulse current"])
            pylab.title('Effective Current ' + "{:6.3f}".format(effectiveCurrent) + u' µA/cm^2')
            pylab.ylabel('Membrane Potential (mV)')
            pylab.xlabel('Time (ms)')
            pylab.ylim([-20,120])
            pylab.xlim([0,60])
            lines += [voltageLine, currentLine]

    if i % 25 == 0:
        # don't wanna print a lot cuz it slows things down
        print "Time: " + str(i*dt) + " ms"
        images.append(tuple(lines))


anim = ArtistAnimation(animationFigure, images, interval = 40, blit = True)
print "Saving animation..."
anim.save("currentComparisons.mp4", extra_args=['-vcodec', 'libx264'])
pylab.show()

# pylab.figure()
# pylab.plot(timeLine, Vm, timeLine, I)
# pylab.title('Hodgkin-Huxley Example')
# pylab.ylabel('Membrane Potential (mV)')
# pylab.xlabel('timeLine (msec)')
# pylab.savefig("model.jpg")

