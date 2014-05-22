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
def createCurrent(distance, stimulusCurrent, t1=5, t2=30):
    effectiveCurrent = stimulusCurrent / (2*np.pi * distance**2) # uA/cm2. the current that reaches the axon.
    I = np.zeros(len(timeLine)) # timecourse of current
    for i, t in enumerate(timeLine):
        if t1 <= t <= t2:
            I[i] = effectiveCurrent
    return I

# Main loop
animationFigure = pylab.figure(figsize=(16, 9))
images = []
for k in range(1, 10):
    pylab.subplot(3, 3, k)

    # create a pulsey current
    distance = k*0.1 # cm. The distance between the axon and the current
    stimulusCurrent = 10. # uA. the current applied at the surface
    I = createCurrent(distance, stimulusCurrent)

    # for the plot's title
    effectiveCurrent = stimulusCurrent / (2*np.pi * distance**2) # uA/cm2. the current that reaches the axon.

    axon = AxonNode(distance) # Create an axon DUN DUN DUNNN

    for i in range(1, len(timeLine)):
        sodiumConductance    = axon.params["gBarNa"] * (axon.m**3) * axon.h
        potassiumConductance = axon.params["gBarK"]  * (axon.n**4)
        leakageConductance   = axon.params["gBarL"]

        # integrate the equations on m, h, and n
        axon.m += (axon.alphaM(axon.Vm[i-1]) * (1 - axon.m) - axon.betaM(axon.Vm[i-1])*axon.m) * dt
        axon.h += (axon.alphaH(axon.Vm[i-1]) * (1 - axon.h) - axon.betaH(axon.Vm[i-1])*axon.h) * dt
        axon.n += (axon.alphaN(axon.Vm[i-1]) * (1 - axon.n) - axon.betaN(axon.Vm[i-1])*axon.n) * dt

        # now integrate the changes in V
        sodiumCurrent = sodiumConductance * (axon.Vm[i-1] - axon.params["sodiumPotential"])
        potassiumCurrent = potassiumConductance * (axon.Vm[i-1] - axon.params["potassiumPotential"])
        leakageCurrent = leakageConductance * (axon.Vm[i-1] - axon.params["leakagePotential"])
        axon.Vm.append(axon.Vm[i-1] + (I[i-1] - sodiumCurrent - potassiumCurrent - leakageCurrent) * dt / axon.params["Cm"])

        # update status
        if i % 25 == 0:
            # update status every 25 frames (about a second of video).
            # don't wanna print a lot cuz it slows things down
            # print "Time: " + str(i*dt) + ", Saving frame " + str(i)

            # plot a frame of the graph
            voltageLine, currentLine = pylab.plot(timeLine[:i+1], axon.Vm[:i+1], 'b-', timeLine[:i+1], I[:i+1], 'g-')
            pylab.legend([voltageLine, currentLine], ["Response from cell", "Impulse current"])
            pylab.title('Effective Current ' + "{:6.3f}".format(effectiveCurrent) + u' µA/cm^2')
            pylab.ylabel('Membrane Potential (mV)')
            pylab.xlabel('Time (ms)')
            pylab.ylim([-20,120])
            pylab.xlim([0,60])
            images.append((voltageLine, currentLine))

anim = ArtistAnimation(animationFigure, images, interval = 50, blit = True)
print "Saving animation..."
pylab.show()
# anim.save("model.mp4", dpi=200, extra_args=['-vcodec', 'libx264'])

# pylab.figure()
# pylab.plot(timeLine, Vm, timeLine, I)
# pylab.title('Hodgkin-Huxley Example')
# pylab.ylabel('Membrane Potential (mV)')
# pylab.xlabel('timeLine (msec)')
# pylab.savefig("model.jpg")

