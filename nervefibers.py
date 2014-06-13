#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pylab
import random

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
            "restingVoltage"     : 0.0,      # V_rest (mv)
            "cm"                 : 1.0,      # µF/cm² membrane capacitance per unit area
            "gBarNa"             : 120.0,    # mS/cm² sodium conductance per unit area
            "gBarK"              : 36.0,     # mS/cm² potassium conductance per unit area
            "gBarL"              : 0.25,     # mS/cm² leakage current conductance per unit area
            "sodiumPotential"    : 115.5,    # mV
            "potassiumPotential" : -11.5,    # mv
            "leakagePotential"   : 11.1,     # mV
            "externalResistivity": 300.0,    # Ω•cm
            "internalResistivity": 110.0     # Ω•cm also called axoplasm resistivity
        }

        params["Cm"] = params["cm"] * np.pi * diameter * length # membrane capacitance (µF)
        params["Ga"] = (np.pi*diameter**2) / (4*params["internalResistivity"] * internodalLength) # axial conductance (S)

        self.Vm    = [params["restingVoltage"]] # The axon node's membrane potential
        self.m     = mInf(params["restingVoltage"])
        self.h     = hInf(params["restingVoltage"])
        self.n     = nInf(params["restingVoltage"])

    # integrate response to stimulus current `stimulus`
    def step(self, stimulus, leftNode, rightNode, dt):
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
        sodiumCurrent    = sodiumConductance    * (self.Vm[lastVm] - self.params["sodiumPotential"])
        potassiumCurrent = potassiumConductance * (self.Vm[lastVm] - self.params["potassiumPotential"])
        leakageCurrent   = leakageConductance   * (self.Vm[lastVm] - self.params["leakagePotential"])

        # MCNEALLLLLL (1976)
        def extV(stimulus, distance): # the external potential
            if distance == 0:
                return 0.0
            else:
                return (self.params["externalResistivity"] * stimulus) / (4 * np.pi * distance)

        neighbourPotential = leftNode["V"] + rightNode["V"] - (2 * self.Vm[lastVm]) # V_n-1 + V_n+1 - 2Vn
        neighbourExtPotential = extV(I, leftNode["d"]) + extV(I, rightNode["d"]) - (2 * extV(I, self.distance))
        ionicCurrent = -np.pi * self.diameter * self.length * (sodiumCurrent + potassiumCurrent + leakageCurrent)

        newV = self.Vm[lastVm]
        newV += (dt / self.params["Cm"]) * (self.params["Ga"] * (neighbourPotential + neighbourExtPotential + ionicCurrent))
        self.Vm.append(newV)

### Simulayshun
class NerveBundleSimulation:
    def __init__(self, T=55, dt=0.025):
        self.T = T
        self.dt = dt
        self.timeLine = np.arange(0, T+dt, dt)

    def simulate(self, nerve, stimulusCurrent):
        for t in range(1, len(self.timeLine)):
            if self.timeLine[t] % 1 == 0.0:
                print "Simulation Time: " + str(self.timeLine[t])

            for i, fiber in enumerate(nerve["fibers"]): # for each nerve fiber
                for k, axonNode in enumerate(fiber.axonNodes): # for each node in the fiber
                    distance = axonNode.distance
                    effectiveCurrent = getCurrent(t*self.dt, stimulusCurrent["magnitude"])

                    lastStep = len(axonNode.Vm) - 1
                    leftNode = {"V": 0, "d": 0}
                    rightNode = {"V": 0, "d": 0}

                    if (k-1) > -1:
                        leftNode["V"] = fiber.axonNodes[k-1].Vm[lastStep]
                        leftNode["d"] = fiber.axonNodes[k-1].distance

                    if (k+1) < len(fiber.axonNodes):
                        rightNode["V"] = fiber.axonNodes[k+1].Vm[lastStep]
                        rightNode["d"] = fiber.axonNodes[k+1].distance

                    # step the current axon forward IN TIIIME ♪♪
                    # print "Stepping axon #" + str(k) + " in fiber #" + str(i)
                    axonNode.step(effectiveCurrent, leftNode, rightNode, self.dt)

# represents a square wave current strimulus
def getCurrent(t, current, tPulseStart=5, pulseWidth=25):
    if tPulseStart <= t <= (tPulseStart + pulseWidth):
        return current
    else:
        return 0.0

def getDistance(x1, y1, z1, x2, y2, z2):
    x = abs(x1 - x2)
    y = abs(y1 - y2)
    z = abs(z1 - z2)
    return np.sqrt(x**2 + y**2 + z**2)

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
                distBetweenFibers = getDistance(x, y, 0, fiber.x, fiber.y, 0)
                if distBetweenFibers < (fiber.diameter/2.0 + diameter/2.0):
                    if distBetweenFibers < fiber.diameter/2.0 or distBetweenFibers < diameter/2.0:
                        # fiber's center is inside another fiber. push it out
                        # this is vector stuff
                        direction = [x - fiber.x, y - fiber.y]
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

        return NerveFiber(x, y, diameter)

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

def plotPositions(plotStimulusPos=True):
    x = []
    y = []
    for i in range(0, len(nerve["fibers"])):
        x.append(nerve["fibers"][i].x)
        y.append(nerve["fibers"][i].y)

    if plotStimulusPos:
        x.append(stimulusCurrent["x"])
        y.append(stimulusCurrent["y"])

    pylab.scatter(x, y, color='r', marker='o', s=0)
    pylab.ylabel('y (cm)')
    pylab.xlabel('x (cm)')
    pylab.axis('equal')
    for i in range(0, len(nerve["fibers"])):
        circle = pylab.Circle((x[i],y[i]), nerve["fibers"][i].diameter/2.0, alpha=0.5)
        pylab.gca().add_artist(circle)

    if plotStimulusPos:
        pylab.axhline(y = stimulusCurrent["y"], color='k', linestyle='--')
        pylab.text(stimulusCurrent["x"]-4, stimulusCurrent["y"]-0.5, "inside")
        pylab.text(stimulusCurrent["x"]-4, stimulusCurrent["y"]+0.2, "outside")
    pylab.show()

# plots the membrane potential of the axons closest to the stimulus
def plotClosestAxons():
    for i in range(0, len(nerve["fibers"])):
        curr = []
        node = nerve["fibers"][i].axonNodes[0]
        for j in range(0, len(simulation.timeLine)):
            curr.append(getCurrent(j*dt, stimulusCurrent["magnitude"]))

        print "plotting axon #" + str(i) + "..."
        pylab.figure()
        pylab.plot(simulation.timeLine, node.Vm, simulation.timeLine, curr)
        pylab.title('Axon #' + str(i) + ": Distance = " + str(node.distance) + " cm")
        pylab.ylabel('Membrane Potential (mV)')
        pylab.xlabel('Time (msec)')
        pylab.savefig("axons/axon" + str(i) + ".jpg")
        pylab.close()

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
    pylab.close();

##############
# Start Script
##############

# Current Stimulus
stimulusCurrent = {
    "magnitude" : 10000, # µA. the current applied at the surface
    "x"         : 0,     # cm
    "y"         : 10,    # cm
    "z"         : 0      # cm
}

# the nerve is a bundle of nerve fibers. Nerve fibers are rods of connected axons.
nerve = {
    "numFibers"   : 50,
    "numNodes"    : 10,    # the number of axon nodes each fiber has
    "fibers"      : [],
    "radius"      : 0.2,   # cm
    "x"           : 0.0,   # cm
    "y"           : 0.0,   # cm
    "z"           : 0.0,   # cm
    "minFiberDiam" : 0.01, # cm
    "maxFiberDiam" : 0.05  # cm
}

class NerveFiber:
    """
    Nerve fibers are myelinated axons that have multiple connected axon nodes (areas of the axon that aren't covered
    by myelin). Axonal Length is in centimetres.
    """
    def __init__(self, x, y, diameter, axonalLength=0.00025, numNodes=10):
        self.x = x
        self.y = y
        self.diameter = diameter
        self.internodalLength = internodalLength = diameter * 100 # McNeal (1976)

        axonalDiameter = self.axonalDiameter = 0.7 * diameter
        axonNodes = self.axonNodes = []
        for i in range(0, numNodes):
            axonNode = AxonPositionNode(i*(axonalLength+internodalLength), axonalDiameter, axonalLength, i, internodalLength)

            nodePos = (self.x, self.y, axonNode.z)  # (x, y, z)
            currPos = (stimulusCurrent["x"], stimulusCurrent["y"], stimulusCurrent["z"]) # (x, y, z)

            # distance from stimulus
            distance = getDistance(nodePos[0], nodePos[1], nodePos[2], currPos[0], currPos[1], currPos[2])
            axonNode.distance = distance
            axonNodes.append(axonNode)

# Create and place the axons
print "Creating bundle of fibers..."
for i in range(0, nerve["numFibers"]):
    succeeded = placeFiberInNerve(nerve)
    if not succeeded:
        break
print "Placed " + str(len(nerve["fibers"])) + " fibers."

# fiber = nerve["fibers"][0]
# for i, axonNode in enumerate(fiber.axonNodes):
#     print "=================" + str(i)
#     print "Node's z position: " + str(axonNode.z)
#     print "Node's diameter: " + str(axonNode.diameter)
#     print "Node's index: " + str(axonNode.index)


T    = 55    # ms
dt   = 0.025 # ms
simulation = NerveBundleSimulation(T, dt)
print "Starting simulation..."
simulation.simulate(nerve, stimulusCurrent) # modifies `nerve`
plotClosestAxons()
print "Done."
