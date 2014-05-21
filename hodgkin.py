#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pylab
import time
from matplotlib.animation import ArtistAnimation

# sauce: http://www.neurdon.com/2011/01/26/neural-modeling-with-python-part-2/

# Potassium (K) Channel
alphaN = np.vectorize(lambda v: 0.01*(10 - v) / (np.exp((10-v)/10) - 1) if v != 10 else 0.1)
betaN  = lambda v: 0.125 * np.exp(-v/80)
nInf   = lambda v: alphaN(v)/(alphaN(v) + betaN(v))

# Sodium (Na) Channel (activating)
alphaM = np.vectorize(lambda v: 0.1*(25-v) / (np.exp((25-v)/10) - 1) if v!= 25 else 1)
betaM = lambda v: 4 * np.exp(-v/18)
mInf = lambda v: alphaM(v)/(alphaM(v) + betaM(v))

# Sodium (Na) Channel (inactivating)
alphaH = lambda v: 0.07 * np.exp(-v/20)
betaH  = lambda v: 1/(np.exp((30-v)/10) + 1)
hInf   = lambda v: alphaH(v)/(alphaH(v) + betaH(v))

# Channel Activity
v = np.arange(-50, 151) # millivolts
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

# Hodgkin-Huxley Parametahs (from the papah!)
restingVoltage     = 0      # V_rest (mv)
Cm                 = 1      # uF/cm2
gBarNa             = 120    # mS/cm2
gBarK              = 36     # mS/cm2
gBarL              = 0.3    # mS/cm2
sodiumPotential    = 115    # mV
potassiumPotential = -12    # mv
leakagePotential   = 10.613 # mV

Vm    = np.zeros(len(timeLine)) # The membrane potential we wanna find
Vm[0] = restingVoltage
m     = mInf(restingVoltage)
h     = hInf(restingVoltage)
n     = nInf(restingVoltage)

# Current Stimulus
def createCurrent(distance, stimulusCurrent, t1=5, t2=30):
    effectiveCurrent = stimulusCurrent / (2*np.pi * distance**2) # uA/cm2. the current that reaches the axon.
    I = np.zeros(len(timeLine)) # timecourse of current
    for i, t in enumerate(timeLine):
        if t1 <= t <= t2:
            I[i] = effectiveCurrent
    return I


# Main loop
animationFigure = pylab.figure()
images = []
for i in range(1, 10):
    pylab.subplot(3, 3, i)

    # create a pulsey current
    distance = i*0.1 # cm. The distance between the axon and the current
    stimulusCurrent = 10. # uA. the current applied at the surface
    I = createCurrent(distance, stimulusCurrent)

    # for the plot's title
    effectiveCurrent = stimulusCurrent / (2*np.pi * distance**2) # uA/cm2. the current that reaches the axon.

    for i in range(1, len(timeLine)):
        sodiumConductance    = gBarNa * (m**3) * h
        potassiumConductance = gBarK  * (n**4)
        leakageConductance   = gBarL

        # integrate the equations on m, h, and n
        m += (alphaM(Vm[i-1]) * (1 - m) - betaM(Vm[i-1])*m) * dt
        h += (alphaH(Vm[i-1]) * (1 - h) - betaH(Vm[i-1])*h) * dt
        n += (alphaN(Vm[i-1]) * (1 - n) - betaN(Vm[i-1])*n) * dt

        # now integrate the changes in V
        sodiumCurrent = sodiumConductance * (Vm[i-1] - sodiumPotential)
        potassiumCurrent = potassiumConductance * (Vm[i-1] - potassiumPotential)
        leakageCurrent = leakageConductance * (Vm[i-1] - leakagePotential)
        Vm[i] = Vm[i-1] + (I[i-1] - sodiumCurrent - potassiumCurrent - leakageCurrent) * dt / Cm

        # update status
        if i % 25 == 0:
            # update status every 25 frames (about a second of video).
            # don't wanna print a lot cuz it slows things down
            # print "Time: " + str(i*dt) + ", Saving frame " + str(i)

            # plot a frame of the graph
            voltageLine, currentLine = pylab.plot(timeLine[:i+1], Vm[:i+1], 'b-', timeLine[:i+1], I[:i+1], 'g-')
            pylab.legend([voltageLine, currentLine], ["Response from cell", "Impulse current"])
            pylab.title('Effective Current ' + "{:6.3f}".format(effectiveCurrent) + u' µA')
            pylab.ylabel('Membrane Potential (mV)')
            pylab.xlabel('Time (ms)')
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

