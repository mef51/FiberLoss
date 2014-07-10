#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pylab
import random
import log

USE_UNITS = False

if USE_UNITS:
    from unum.units import *
    from unum import Unum
else:
    from fastunits import *

## Setup Units
mV = Unum.unit('mV', 10**-3 * V) # millivolts
mF = Unum.unit('mF', 10**-3 * F) # milliFarads
uF = Unum.unit('uF', 10**-6 * F) # microFarads

mS = Unum.unit('mS', 10**-3 * S) # milliSiemens
mohm = Unum.unit('mohm', 10**-3 * ohm) # milliohms

# strips the unit from a number and returns just the number. `unit` is a Unum unit.
def mag(v, unit):
    return float(v/unit)

class AxonPositionNode:
    """
    A node of Ranvier on an Axon, as modelled by Hodgkin and Huxley in 1952 .
    This class is meant for creating axons with a specific location
    """
    def __init__(self, z, diameter, length, index, internodalLength):
        # position
        self.z = z               # position down the fiber
        self.diameter = diameter # the diameter of the node
        self.length = length     # the length of the node

        # each axon in a node is labelled with a number (n). the axon closest to the stimulus is numbered n = 0
        self.index = index

        # Hodgkin-Huxley Parametahs (from the papah!)
        params = self.params = {
            "restingVoltage"     : 0       *mV,         # V_rest (mv)
            "cm"                 : 1       *mF/(cm**2), # mF/cm² membrane capacitance per unit area
            "gBarNa"             : 120.0   *mS/(cm**2), # mS/cm² sodium conductance per unit area
            "gBarK"              : 36.0    *mS/(cm**2), # mS/cm² potassium conductance per unit area
            "gBarL"              : 0.25    *mS/(cm**2), # mS/cm² leakage current conductance per unit area
            "sodiumPotential"    : (50.5+70)    *mV,         # mV
            "potassiumPotential" : (-77.0+70)   *mV,         # mv
            "leakagePotential"   : (-54.4+70)   *mV,         # mV
            "externalResistivity": 300.0e3 *mohm*cm,    # mΩ•cm
            "internalResistivity": 110.0e3 *mohm*cm     # mΩ•cm also called axoplasm resistivity
        }

        ###### Potassium (K) Channel
        def alphaN(v):
            if v == 35*mV: return 0.2 * (1/ms)
            a = 1 - np.exp(float((35*mV - v)/(10.0*mV)))
            a = a**-1
            a *= 0.02*(1/(mV*ms))*(v-35*mV)
            return a
        self.alphaN = alphaN = np.vectorize(alphaN)

        def betaN(v):
            if v == 10*mV: return 0.5 * (1/ms)
            b = (1 - np.exp(float((v - 10*mV)/(10.0*mV)))) ** -1
            b *= 0.05*(1/(mV*ms))*(10*mV - v)
            return b
        self.betaN = betaN = np.vectorize(betaN)
        nInf = self.nInf   = lambda v: alphaN(v)/(alphaN(v) + betaN(v))

        ###### Sodium (Na) Channel (activating)
        def alphaM(v):
            if v == 22*mV: return 1.08 * (1/ms)
            a = (1 - np.exp(float((22*mV - v)/(3.0*mV)))) ** -1
            a *= 0.36*(1/(mV*ms))*(v - 22*mV)
            return a
        self.alphaM = alphaM = np.vectorize(alphaM)

        def betaM(v):
            if v == 13*mV: return 8.0 * (1/ms)
            b = (1 - np.exp(float((v - 13*mV)/(20.0*mV)))) ** -1
            b *= 0.4*(1/(mV*ms))*(13*mV - v)
            return b
        self.betaM = betaM = np.vectorize(betaM)
        mInf = self.mInf = lambda v: alphaM(v)/(alphaM(v) + betaM(v))

        ###### Sodium (Na) Channel (inactivating)
        def alphaH(v):
            if v == -10*mV: return 0.6 * (1/ms)
            a = (1 - np.exp(float((v + 10*mV)/(6.0*mV)))) ** -1
            a *= 0.1 * (1/(mV*ms)) * (-10*mV - v)
            return a
        self.alphaH = alphaH = np.vectorize(alphaH)

        def betaH(v):
            b = (1 + np.exp(float((45*mV - v)/(10.0*mV)))) ** -1
            b *= 4.5 * (1/ms)
            return b
        betaH = self.betaH = np.vectorize(betaH)
        hInf = self.hInf   = lambda v: alphaH(v)/(alphaH(v) + betaH(v))

        def sodiumCurrent(V, m, h):
            return self.params["gBarNa"] * (m**3) * h * (V - self.params["sodiumPotential"])
        self.sodiumCurrent = sodiumCurrent

        def potassiumCurrent(V, n):
            return self.params["gBarK"]  * (n**4) * (V - self.params["potassiumPotential"])
        self.potassiumCurrent = potassiumCurrent

        def leakageCurrent(V):
            return self.params["gBarL"]  * (V - self.params["leakagePotential"])
        self.leakageCurrent = leakageCurrent

        # MCNEALLLLLL (1976)
        def extV(stimulus, distance): # the external potential
            if mag(distance, cm) == 0:
                return 0.0*mV
            else:
                V = (self.params["externalResistivity"] * stimulus) / (4 * np.pi * distance)
                return V
        self.extV = extV

        log.infoVar(diameter, 'diameter')
        log.infoVar(length, 'length')
        log.infoVar(np.pi, 'pi')
        log.infoVar(params["cm"], 'cm')

        params["Cm"] = params["cm"] * np.pi * diameter * length # membrane capacitance (mF)
        params["Ga"] = (np.pi*diameter**2) / (4*params["internalResistivity"] * internodalLength) # axial conductance (mS)

        self.Vm = [params["restingVoltage"]] # The axon node's membrane potential
        self.m  = [mInf(params["restingVoltage"])]
        self.h  = [hInf(params["restingVoltage"])]
        self.n  = [nInf(params["restingVoltage"])]

    # integrate response to stimulus current `stimulus`
    def step(self, stimulus, leftNode, rightNode, dt):
        I = stimulus # I[i-1]
        extV = self.extV

        def last(l):
            return l[len(l) - 1]

        V, m, h, n = last(self.Vm), last(self.m), last(self.h), last(self.n)

        iNa = self.sodiumCurrent(V, m, h)
        iK  = self.potassiumCurrent(V, n)
        iL  = self.leakageCurrent(V)

        log.infoVar(m, "self.m")
        log.infoVar(h, "self.h")
        log.infoVar(n, "self.n")

        def checkGates(gate):
            if not (0 < gate < 1):
                log.error("Gating variable is out of range! D: " + str(gate))
                exit()

        for gate in [m, h, n]:
            checkGates(gate)

        # integrate the equations on m, h, and n
        newM = (self.alphaM(V) * (1 - m) - self.betaM(V)*m) * dt + m
        newH = (self.alphaH(V) * (1 - h) - self.betaH(V)*h) * dt + h
        newN = (self.alphaN(V) * (1 - n) - self.betaN(V)*n) * dt + n

        neighbourPotential = leftNode["V"] + rightNode["V"] - (2 * V) # V_n-1 + V_n+1 - 2Vn
        neighbourExtPotential = extV(I, leftNode["d"]) + extV(I, rightNode["d"]) - (2 * extV(I, self.distance))

        surroundingCurrent = self.params["Ga"] * (neighbourPotential + neighbourExtPotential)
        ionicCurrent = np.pi * self.diameter * self.length * (iNa + iK + iL)

        log.infoVar(neighbourPotential, "neighbourPotential")
        log.infoVar(neighbourExtPotential, "neighbourExtPotential")
        log.infoVar(iNa, "sodiumCurrent")
        log.infoVar(iK, "potassiumCurrent")
        log.infoVar(iL, "leakageCurrent")
        log.infoVar(V, "lastVm")

        # now integrate the changes in V
        newV = (dt / self.params["Cm"]) * (surroundingCurrent - ionicCurrent) + V

        log.infoVar(self.params["Cm"], "Cm")
        log.infoVar(surroundingCurrent, "surroundingCurrent")
        log.infoVar(ionicCurrent, "ionicCurrent")
        log.infoVar(newV, "newV")

        self.m.append(newM)
        self.h.append(newH)
        self.n.append(newN)
        self.Vm.append(newV)

    def plotAlphaBetaFunctions(self):
        v = np.arange(-75, 125) # millivolts
        pylab.figure()
        pylab.xlim([-75, 125])
        pylab.plot(v, self.alphaM(v), v, self.alphaH(v), v, self.alphaN(v), v, self.betaM(v), v, self.betaH(v), v, self.betaN(v))
        pylab.legend(('alphaM', 'alphaH', 'alphaN', 'betaM', 'betaH', 'betaN'))
        pylab.title('Alpha and Beta Functions')
        pylab.ylabel(u'Rate Constant (ms^-1)')
        pylab.xlabel('Voltage (mV)')
        pylab.savefig('alphaBetaFunctions.jpg')

    def plotCurrentsVoltagesAndGates(self, timeLine, stimulusCurrent, fiberNum, plotStimulus=True):
        vSol, mSol, hSol, nSol = self.Vm, self.m, self.h, self.n
        extPotentialSol = [self.extV(stimulusCurrent[i], self.distance) for i, t in enumerate(timeLine)]

        # current solutions
        iNaSol = [self.sodiumCurrent(vSol[i], mSol[i], hSol[i]) for i in range(0, len(timeLine))]
        iKSol  = [self.potassiumCurrent(vSol[i], nSol[i]) for i in range(0, len(timeLine))]
        iLSol  = [self.leakageCurrent(vSol[i]) for i in range(0, len(timeLine))]

        # strip out units
        vSol   = [mag(val, mV) for val in vSol]
        mSol   = [mag(val, 1) for val in mSol]
        hSol   = [mag(val, 1) for val in hSol]
        nSol   = [mag(val, 1) for val in nSol]
        iNaSol = [mag(val, mA/(cm*cm)) for val in iNaSol]
        iKSol  = [mag(val, mA/(cm*cm)) for val in iKSol]
        iLSol  = [mag(val, mA/(cm*cm)) for val in iLSol]
        extPotentialSol = [mag(val, mV) for val in extPotentialSol]

        pylab.subplot(1, 2, 1)
        if plotStimulus:
            pylab.plot(timeLine, vSol, timeLine, extPotentialSol)
        else:
            pylab.plot(timeLine, vSol)

        d = mag(self.distance, cm)
        pylab.title("Membrane Voltage of node #" + str(self.index) + ": d = " + str("{0:.2f}".format(d)) + "cm")
        pylab.legend(('V', 'Stimulus'))
        pylab.ylabel('V (mV)')
        pylab.xlabel('Time (ms)')
        pylab.grid()

        pylab.subplot(2, 2, 2)
        pylab.plot(timeLine, iNaSol, timeLine, iKSol, timeLine, iLSol)
        pylab.title("Ionic Currents")
        pylab.legend(('iNa', 'iK', 'iL'))
        pylab.ylabel('Current (mA)')
        pylab.xlabel('Time (ms)')
        pylab.grid()

        pylab.subplot(2, 2, 4)
        pylab.plot(timeLine, mSol, timeLine, hSol, timeLine, nSol)
        pylab.title("m, h, and n")
        pylab.legend(('m', 'h', 'n'))
        pylab.ylabel('Probability')
        pylab.xlabel('Time (ms)')
        pylab.grid()

        pylab.tight_layout()
        pylab.savefig("graphs/axons/axon" + str(self.index) + "fiber" + str(fiberNum) + ".jpg")
        pylab.close()

class NerveFiber:
    """
    Nerve fibers are myelinated axons that have multiple connected axon nodes (areas of the axon that aren't covered
    by myelin). Axonal Length is in centimetres.
    """
    def __init__(self, x, y, diameter, numNodes, axonalLength):
        self.x = x
        self.y = y
        self.diameter = diameter
        self.internodalLength = internodalLength = diameter * 100 # McNeal (1976)

        axonalDiameter = self.axonalDiameter = 0.7 * diameter
        axonNodes = self.axonNodes = []
        for i in range(0, numNodes):
            index = i - int(numNodes/2)
            z = index*(axonalLength+internodalLength)
            axonNode = AxonPositionNode(z, axonalDiameter, axonalLength, index, internodalLength)

            nodePos = (self.x, self.y, axonNode.z)  # (x, y, z)
            currPos = (stimulusCurrent["x"], stimulusCurrent["y"], stimulusCurrent["z"]) # (x, y, z)

            # distance from stimulus
            distance = getDistance(nodePos[0], nodePos[1], nodePos[2], currPos[0], currPos[1], currPos[2])
            axonNode.distance = distance
            axonNodes.append(axonNode)

### Simulayshun
class NerveBundleSimulation:
    def __init__(self, T=55*ms, dt=0.025*ms):
        self.T = T
        self.dt = dt
        self.timeLine = np.arange(0, mag((T+dt), ms), mag(dt, ms))

    def simulate(self, nerve, stimulusCurrent):
        for t in range(1, len(self.timeLine)):
            if self.timeLine[t] % 1 == 0.0:
                print "Simulation Time: " + str(self.timeLine[t])

            for i, fiber in enumerate(nerve["fibers"]): # for each nerve fiber
                for k, axonNode in enumerate(fiber.axonNodes): # for each node in the fiber
                    distance = axonNode.distance
                    effectiveCurrent = getCurrent(t*self.dt, stimulusCurrent["magnitude"])

                    lastStep = len(axonNode.Vm) - 1
                    leftNode  = {"V": 0.0*mV, "d": 0.0*cm}
                    rightNode = {"V": 0.0*mV, "d": 0.0*cm}

                    if (k-1) > -1:
                        leftNode["V"] = fiber.axonNodes[k-1].Vm[lastStep]
                        leftNode["d"] = fiber.axonNodes[k-1].distance

                    if (k+1) < len(fiber.axonNodes):
                        rightNode["V"] = fiber.axonNodes[k+1].Vm[lastStep]
                        rightNode["d"] = fiber.axonNodes[k+1].distance

                    # step the current axon forward IN TIIIME ♪♪
                    # print "Stepping axon #" + str(k) + " in fiber #" + str(i)
                    log.info("========")
                    log.infoVar(t*dt, 'time')
                    log.infoVar(axonNode.index, 'node')
                    log.infoVar(leftNode["V"], "leftVoltage")
                    log.infoVar(rightNode["V"], "rightVoltage")
                    axonNode.step(effectiveCurrent, leftNode, rightNode, self.dt)

# represents a square wave current strimulus
def getCurrent(t, current, tPulseStart=0*ms, pulseWidth=25*ms):
    if tPulseStart <= t <= (tPulseStart + pulseWidth):
        return current
    else:
        return 0.0*mA

def getDistance(x1, y1, z1, x2, y2, z2):
    x = abs(x1 - x2)
    y = abs(y1 - y2)
    z = abs(z1 - z2)
    return np.sqrt(mag((x**2 + y**2 + z**2), (cm*cm))) * (cm)

# Places a nerve fiber and makes sure it doesn't overlap with any other nerve fibers in the nerve bundle
# returns true if it succeeds, false otherwise
def placeFiberInNerve(nerve, maxAttempts = 1000):
    def placeFiber(maxAttempts = 1000):
        x = random.uniform(-nerve["radius"], nerve["radius"]) + nerve["x"]
        y = random.uniform(-nerve["radius"], nerve["radius"]) + nerve["y"]
        diameter = random.uniform(nerve["minFiberDiam"], nerve["maxFiberDiam"])

        # make sure axon is in the nerve
        while (x**2 + y**2) > nerve["radius"]**2:
            x = random.uniform(-nerve["radius"], nerve["radius"])
            y = random.uniform(-nerve["radius"], nerve["radius"])

        # push around and shrink the new axon until it fits
        recheck = True
        numRechecks = 0
        while recheck:
            if numRechecks > maxAttempts:
                print "Couldn't place " + str(nerve["numFibers"]) + " fibers after " + str(maxAttempts) + " attempts!"
                print "Giving up at " + str(len(nerve["fibers"])) + " fibers."
                return None

            recheck = False
            for k, fiber in enumerate(nerve["fibers"]):
                distBetweenFibers = getDistance(x, y, 0*cm, fiber.x, fiber.y, 0*cm)
                if distBetweenFibers < (fiber.diameter/2.0 + diameter/2.0):
                    if distBetweenFibers < fiber.diameter/2.0 or distBetweenFibers < diameter/2.0:
                        # fiber's center is inside another fiber. push it out
                        # this is vector stuff
                        direction = [mag(x - fiber.x, cm), mag(y - fiber.y, cm)]
                        dirMag = np.sqrt(direction[0]**2 + direction[1]**2)
                        # normalize this vector so that we can use it as a direction
                        direction[0] = direction[0] /dirMag
                        direction[1] = direction[1] /dirMag

                        x += direction[0] * (fiber.diameter/2.0 - distBetweenFibers + diameter/2.0)
                        y += direction[1] * (fiber.diameter/2.0 - distBetweenFibers + diameter/2.0)
                        recheck = True
                        numRechecks += 1
                        break
                    else:
                        # fiber is too big
                        amountToShrink = (fiber.diameter/2.0 + diameter/2.0) - distBetweenFibers
                        diameter -= amountToShrink * 2

        return NerveFiber(x, y, diameter, nerve["numNodes"], nerve["axonalLength"])

    fiber = placeFiber(maxAttempts)
    # make sure fiber isn't too big or small
    while fiber.diameter < nerve["minFiberDiam"] or fiber.diameter > nerve["maxFiberDiam"]:
        result = placeFiber(maxAttempts)
        if result is None:
            return False
        else:
            fiber = result

    nerve["fibers"].append(fiber)
    return True

##########################
# Some Plotting Functions
##########################

def plotNodePositions(plotIndices=True):
    pylab.figure()
    pylab.title("Node Positions")
    pylab.ylabel('y (cm)')
    pylab.xlabel('z (cm)')
    xSpan = mag(nerve["fibers"][0].internodalLength*nerve["numNodes"], cm)
    ySpan = mag(nerve["fibers"][0].diameter, cm) + abs(mag(nerve["fibers"][0].y, cm))
    pylab.xlim([-xSpan, xSpan])
    pylab.ylim([-ySpan, ySpan])

    for i, fiber in enumerate(nerve["fibers"]):
        for j, node in enumerate(fiber.axonNodes):
            newZ = mag(node.z, cm)
            newY = mag(fiber.y, cm)

            rectangle = pylab.Rectangle((newZ,newY), mag(node.length, cm), mag(node.diameter, cm), alpha=0.5)
            pylab.gca().add_artist(rectangle)
            if plotIndices:
                pylab.text(newZ, newY + mag(node.diameter, cm), str(node.index))

    pylab.show()

def plotCrossSectionPositions(plotStimulusPos=True):
    x = []
    y = []
    for i in range(0, len(nerve["fibers"])):
        x.append(mag(nerve["fibers"][i].x, cm))
        y.append(mag(nerve["fibers"][i].y, cm))

    if plotStimulusPos:
        x.append(mag(stimulusCurrent["x"], cm))
        y.append(mag(stimulusCurrent["y"], cm))

    pylab.scatter(x, y, color='r', marker='o', s=0)
    pylab.ylabel('y (cm)')
    pylab.xlabel('x (cm)')
    pylab.axis('equal')
    for i in range(0, len(nerve["fibers"])):
        circle = pylab.Circle((x[i],y[i]), mag(nerve["fibers"][i].diameter/(2.0), cm), alpha=0.5)
        pylab.gca().add_artist(circle)

    if plotStimulusPos:
        pylab.axhline(y = mag(stimulusCurrent["y"], cm), color='k', linestyle='--')
        pylab.text(mag(stimulusCurrent["x"], cm)-4, mag(stimulusCurrent["y"], cm)-0.5, "inside")
        pylab.text(mag(stimulusCurrent["x"], cm)-4, mag(stimulusCurrent["y"], cm)+0.2, "outside")
    pylab.show()

def plotMembranePotential(current, node, fiberNum, showFigure=False):
    pylab.figure()

    # strip out units
    Vm = [mag(v, mV) for v in node.Vm]
    curr = [mag(c, mA) for c in current]

    # plot
    pylab.plot(simulation.timeLine, Vm, simulation.timeLine, curr)
    pylab.title('Axon #' + str(node.index) + ": Distance = " + str(node.distance) + " cm")
    pylab.ylabel('Membrane Potential (mV)')
    pylab.xlabel('Time (msec)')
    pylab.savefig("graphs/axons/axon" + str(node.index) + "fiber" + str(fiberNum) + ".jpg")
    if showFigure:
        pylab.show()
    pylab.close()

# plots the membrane potential of the axons closest to the stimulus
def plotClosestAxons():
    curr = []
    for j in range(0, len(simulation.timeLine)):
        curr.append(getCurrent(j*dt, stimulusCurrent["magnitude"]))

    for i in range(0, len(nerve["fibers"])):
        node = nerve["fibers"][i].axonNodes[0]
        print "plotting axon #" + str(node.index) + " in fiber #" + str(i) + "..."
        plotMembranePotential(curr, node, i)

def plotMembranePotentialOfNodes(nerve):
    curr = []
    for j in range(0, len(simulation.timeLine)):
        curr.append(getCurrent(j*dt, stimulusCurrent["magnitude"]))

    for i, fiber in enumerate(nerve["fibers"]):
        for k, node in enumerate(fiber.axonNodes):
            print "plotting axon #" + str(node.index) + " in fiber #" + str(i) + "..."
            plotMembranePotential(curr, node, i)

def plotInfoOfNodes(nerve, plotStimulus=True):
    curr = []
    for j in range(0, len(simulation.timeLine)):
        curr.append(getCurrent(j*dt, stimulusCurrent["magnitude"]))

    for fiberNum, fiber in enumerate(nerve["fibers"]):
        for k, node in enumerate(fiber.axonNodes):
            print "plotting axon #" + str(node.index) + " in fiber #" + str(fiberNum) + "..."
            node.plotCurrentsVoltagesAndGates(simulation.timeLine, curr, fiberNum, plotStimulus)

def plotCompoundPotential():
    compoundPotential = []
    for i in range(0, int(T/dt) + 1):
        compoundPotential += [0]

    for i in range(0, len(nerve['fibers'])):
        axon = nerve['fibers'][i]
        for k, v in enumerate(axon.Vm):
            compoundPotential[k] += v

    pylab.figure()
    pylab.plot(simulation.timeLine, compoundPotential)
    pylab.xlabel('Time (ms)')
    pylab.ylabel('Voltage (mV)')
    pylab.title('Sum of all Action Potentials')
    pylab.savefig("sumOfPotentials.jpg")
    pylab.close()

##############
# Start Script
##############

log.logLevel = log.ERROR
# log.logLevel = log.INFO

# Current Stimulus
stimulusCurrent = {
    "magnitude" : -0.3 *mA,    # mA. the current applied at the surface
    "x"         : 0   *cm,    # cm
    "y"         : 0.1 *cm,    # cm
    "z"         : 0   *cm     # cm
}

# the nerve is a bundle of nerve fibers. Nerve fibers are rods of connected axons.
nerve = {
    "numFibers"    : 1,
    "numNodes"     : 1,    # the number of axon nodes each fiber has
    "fibers"       : [],
    "radius"       : 0.2    *cm, # cm
    "x"            : 0.0    *cm, # cm
    "y"            : 0.0    *cm, # cm
    "z"            : 0.0    *cm, # cm
    "minFiberDiam" : 0.0019 *cm, # cm
    "maxFiberDiam" : 0.0021 *cm, # cm
    "axonalLength" : 2.5e-4 *cm # cm
}

# Create and place the axons
print "Creating bundle of fibers..."
for i in range(0, nerve["numFibers"]):
    succeeded = placeFiberInNerve(nerve)
    if not succeeded:
        break
print "Placed " + str(len(nerve["fibers"])) + " fibers."

# plotNodePositions()

T    = 55*ms    # ms
dt   = 0.025*ms # ms
simulation = NerveBundleSimulation(T, dt)
print "Starting simulation..."
simulation.simulate(nerve, stimulusCurrent) # modifies `nerve`
plotInfoOfNodes(nerve, plotStimulus=False)

print "Done."
