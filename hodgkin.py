#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pylab
from matplotlib.animation import ArtistAnimation
from neurons import AxonNode

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
        axon.step(I, dt)

        # update status every 25 frames (about a second of video).
        if i % 25 == 0:
            # get a profile of the current up til now
            curr = []
            for j in range(0, len(timeLine[:i+1])):
                curr.append(getCurrent(j*dt, axon.distance, stimulusCurrent))

            # plot a frame of the graph
            voltageLine, currentLine = pylab.plot(timeLine[:i+1], axon.Vm[:i+1], 'b-', timeLine[:i+1], curr[:i+1], 'g-')
            if k == 1: # have only 1 legend
                pylab.legend([voltageLine, currentLine], ["Response from cell", "Impulse current"])
            pylab.title('Effective Current ' + "{:6.3f}".format(effectiveCurrent) + u' µA/cm$^2$')
            pylab.ylabel('Membrane Potential (mV)')
            pylab.xlabel('Time (ms)')
            pylab.ylim([-20,120])
            pylab.xlim([0,60])
            lines += [voltageLine, currentLine]

    if i % 25 == 0:
        # don't wanna print a lot cuz it slows things down
        print "Time: " + str(i*dt) + " ms"
        images.append(tuple(lines))

pylab.tight_layout()
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

