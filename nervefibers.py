#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pylab
import random
import log
import argparse
import os

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
kohm = Unum.unit('kohm', 10**3 * ohm) # kiloohms

# strips the unit from a number and returns just the number. `unit` is a Unum unit.
def mag(v, unit):
    return float(v/unit)

# little utility for getting the last element of a list
def last(l):
    return l[len(l) - 1]

class AxonPositionNode:
    """
    A node of Ranvier on an Axon, as modelled by Hodgkin and Huxley in 1952 .
    This class is meant for creating axons with a specific location
    """
    def __init__(self, z, diameter, length, index, internodalLength, damagedChannels):
        # position
        self.z = z               # position down the fiber
        self.diameter = diameter # the diameter of the node
        self.length = length     # the length of the node

        # the proportion of sodium ion channels in the membrane that are damaged
        self.damagedChannels = damagedChannels

        # each axon in a node is labelled with a number (n). the axon closest to the stimulus is numbered n = 0
        self.index = index

        # Hodgkin-Huxley Parametahs (from the papah!)
        params = self.params = {
            "restingVoltage"     : -65.5   *mV,         # V_rest (mv)
            "cm"                 : 1       *uF/(cm**2), # mF/cm² membrane capacitance per unit area
            "gBarNa"             : 120.0   *mS/(cm**2), # mS/cm² sodium conductance per unit area
            "gBarK"              : 36.0    *mS/(cm**2), # mS/cm² potassium conductance per unit area
            "gBarL"              : 0.25    *mS/(cm**2), # mS/cm² leakage current conductance per unit area
            "sodiumPotential"    : 50      *mV,         # mV
            "potassiumPotential" : -77     *mV,         # mv
            "leakagePotential"   : -54.4   *mV,         # mV
            "externalResistivity": 300.0   *ohm*cm,     # Ω•cm
            "internalResistivity": 110/3.4 *ohm*cm,     # Ω•cm also called axoplasm resistivity
            "leftShift"          : 35      *mV
        }

        ###### Potassium (K) Channel
        def alphaN(v):
            if v == -55*mV: return (0.1*(1/ms))
            a = -(v + 55*mV)/(10.0*mV)
            a = 1.0/(1 - np.exp(float(a)))
            a *= v + 55*mV
            return a * 0.01*(1/(ms*mV))
        self.alphaN = alphaN = np.vectorize(alphaN)

        def betaN(v):
            b = -(v+65*mV)/(80.0*mV)
            return 0.125 * (1/ms) * np.exp(float(b))
        self.betaN = betaN = np.vectorize(betaN)
        nInf = self.nInf   = lambda v: alphaN(v)/(alphaN(v) + betaN(v))

        ###### Sodium (Na) Channel (activating)
        def alphaM(v):
            if v == -40*mV: return 1.0 * (1/ms)
            a = -(v + 40*mV)/(10.0*mV)
            a = 1.0/(1 - np.exp(float(a)))
            a *= v + 40*mV
            return a * 0.1*(1/(ms*mV))
        self.alphaM = alphaM = np.vectorize(alphaM)

        def betaM(v):
            b = -(v + 65*mV)/(18.0*mV)
            b = 4 * (1/ms) * np.exp(float(b))
            return b
        self.betaM = betaM = np.vectorize(betaM)
        mInf = self.mInf = lambda v: alphaM(v)/(alphaM(v) + betaM(v))

        ###### Sodium (Na) Channel (inactivating)
        def alphaH(v):
            a = -(v + 65*mV)/(20*mV)
            return 0.07 * (1/ms) * np.exp(float(a))
        self.alphaH = alphaH = np.vectorize(alphaH)

        def betaH(v):
            b = -(v + 35*mV)/(10.0*mV)
            return 1.0 * (1/ms) / (1 + np.exp(float(b)))
        betaH = self.betaH = np.vectorize(betaH)
        hInf = self.hInf   = lambda v: alphaH(v)/(alphaH(v) + betaH(v))

        def sodiumCurrent(V, m, h, mLS, hLS):
            current = self.params["gBarNa"] * (V - self.params["sodiumPotential"])
            return  current * (m**3 * h * (1 - self.damagedChannels) + mLS**3 * hLS * self.damagedChannels)
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

        params["Cm"] = params["cm"] * np.pi * diameter * length # membrane capacitance (uF)
        params["Ga"] = (np.pi*diameter**2) / (4*params["internalResistivity"] * internodalLength) # axial conductance (mS)
        params["Kappa"] = 1 / (params["internalResistivity"] * internodalLength) # axial conductance (mS)

        self.Vm = [params["restingVoltage"]] # The axon node's membrane potential
        self.m  = [mInf(params["restingVoltage"])]
        self.h  = [hInf(params["restingVoltage"])]
        self.n  = [nInf(params["restingVoltage"])]

        # Nav-CLS (Boucher 2012), the m and h variables are left shifted for damaged channels
        self.mLS = [mInf(params["restingVoltage"] + self.params["leftShift"])]
        self.hLS = [hInf(params["restingVoltage"] + self.params["leftShift"])]

    # integrate response to stimulus current `stimulus`
    def step(self, stimulus, leftNode, rightNode, dt, exciteCenterOnly=False):
        I = stimulus # I[i-1]
        extV = self.extV
        leftShift = self.params["leftShift"]

        V, m, h, n, mLS, hLS = last(self.Vm), last(self.m), last(self.h), last(self.n), last(self.mLS), last(self.hLS)

        iNa = self.sodiumCurrent(V, m, h, mLS, hLS)
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
        newMLS = (self.alphaM(V + leftShift) * (1 - mLS) - self.betaM(V + leftShift)*mLS) * dt + mLS
        newHLS = (self.alphaH(V + leftShift) * (1 - hLS) - self.betaH(V + leftShift)*hLS) * dt + hLS

        neighbourPotential = leftNode["V"] + rightNode["V"] - (2 * V) # V_n-1 + V_n+1 - 2Vn
        neighbourExtPotential = 0

        if exciteCenterOnly: # only add up the terms for the node with index 0
            if self.index == 0:
                neighbourExtPotential += - (2 * extV(I, self.distance))
            elif leftNode["n"] == 0:
                neighbourExtPotential += extV(I, leftNode["d"])
            elif rightNode["n"] == 0:
                neighbourExtPotential += extV(I, rightNode["d"])
        else:
            neighbourExtPotential += extV(I, leftNode["d"]) + extV(I, rightNode["d"]) - (2 * extV(I, self.distance))

        surroundingCurrent = self.params["Ga"] * (neighbourPotential + neighbourExtPotential)
        ionicCurrent = np.pi * self.diameter * self.length * (iNa + iK + iL)

        log.infoVar(neighbourPotential, "neighbourPotential")
        log.infoVar(neighbourExtPotential, "neighbourExtPotential")
        log.infoVar(iNa, "sodiumCurrent")
        log.infoVar(iK, "potassiumCurrent")
        log.infoVar(iL, "leakageCurrent")
        log.infoVar(V, "lastVm")

        # if self.index != 0: I = 0
        # now integrate the changes in V
        newV = (dt / self.params["Cm"]) * (surroundingCurrent + 0*I - ionicCurrent) + V

        log.infoVar(self.params["Cm"], "Cm")
        log.infoVar(self.params["Ga"], "Ga")
        log.infoVar(self.params["Kappa"], "Kappa")
        log.infoVar(surroundingCurrent, "surroundingCurrent")
        log.infoVar(ionicCurrent, "ionicCurrent")
        log.infoVar(newV, "newV")

        self.m.append(newM)
        self.h.append(newH)
        self.n.append(newN)
        self.mLS.append(newMLS)
        self.hLS.append(newHLS)
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
        saveFigure('alphaBetaFunctions.jpg')

    def plotSteadyStateActivations(self):
        v = np.arange(-150, 50) # millivolts
        leftShift = self.params["leftShift"]
        mInf, hInf, nInf = self.mInf, self.hInf, self.nInf
        pylab.figure()
        pylab.plot(v, mInf(v), v, hInf(v), v, nInf(v), v, mInf(v + leftShift), v, hInf(v + leftShift))
        pylab.legend(('m', 'h', 'n', 'mLS', 'hLS'))
        pylab.title('Steady state values of ion channel gating variables')
        pylab.ylabel('Magnitude')
        pylab.xlabel('Voltage (mV)')
        saveFigure("mhn.jpg")

    def plotCurrentsVoltagesAndGates(self, timeLine, stimulusCurrent, fiberNum, plotStimulus=True, exciteCenterOnly=False):
        vSol, mSol, hSol, nSol, mLSSol, hLSSol = self.Vm, self.m, self.h, self.n, self.mLS, self.hLS

        n = 1
        if exciteCenterOnly and self.index != 0: n = 0
        extPotentialSol = [self.extV(n*stimulusCurrent[i], n*self.distance) for i, t in enumerate(timeLine)]

        # current solutions
        iNaSol = [self.sodiumCurrent(vSol[i], mSol[i], hSol[i], mLSSol[i], hLSSol[i]) for i in range(0, len(timeLine))]
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
        pylab.ylim(-200, 150)
        if plotStimulus:
            pylab.plot(timeLine, vSol, timeLine, extPotentialSol)
        else:
            pylab.plot(timeLine, vSol)

        d = mag(self.distance, cm)
        titleStr = "Vm of node #" + str(self.index) + ": d = " + str("{0:.2f}".format(d)) + "cm. "
        titleStr += "AC = " + str(self.damagedChannels)
        pylab.title(titleStr)
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
        if self.damagedChannels == 0:
            pylab.plot(timeLine, mSol, timeLine, hSol, timeLine, nSol)
            pylab.title("m, h, and n")
            pylab.legend(('m', 'h', 'n'))
        else:
            pylab.plot(timeLine, mSol, timeLine, hSol, timeLine, nSol, timeLine, mLSSol, timeLine, hLSSol)
            pylab.title("m, h, and n with left-shifts")
            pylab.legend(('m', 'h', 'n', 'mLS', 'hLS'))
        pylab.ylabel('Probability')
        pylab.xlabel('Time (ms)')
        pylab.grid()

        pylab.tight_layout()
        saveFigure("graphs/axons/axon" + str(self.index) + "fiber" + str(fiberNum) + ".jpg")
        pylab.close()

class NerveFiber:
    """
    Nerve fibers are myelinated axons that have multiple connected axon nodes (areas of the axon that aren't covered
    by myelin). Axonal Length is in centimetres.
    """
    def __init__(self, x, y, diameter, numNodes, axonalLength, damageMap, stimulusCurrent):
        self.x = x
        self.y = y
        self.diameter = diameter
        self.axonalLength = axonalLength
        self.internodalLength = internodalLength = diameter * 100 # McNeal (1976)
        self.damageMap = damageMap

        axonalDiameter = self.axonalDiameter = 0.7 * diameter
        axonNodes = self.axonNodes = [] # da CANONICAL list of nodes in da fiber bro
        for i in range(0, numNodes):
            index = i - int(numNodes/2)
            z = index*(axonalLength+internodalLength)

            damage = 0.0
            if index in damageMap.keys():
                damage = damageMap[index]

            axonNode = AxonPositionNode(z, axonalDiameter, axonalLength, index, internodalLength, damage)

            nodePos = (self.x, self.y, axonNode.z)  # (x, y, z)
            currPos = (stimulusCurrent["x"], stimulusCurrent["y"], stimulusCurrent["z"]) # (x, y, z)

            # distance from stimulus
            distance = getDistance(nodePos[0], nodePos[1], nodePos[2], currPos[0], currPos[1], currPos[2])
            axonNode.distance = distance
            axonNodes.append(axonNode)

### Simulayshun
class NerveBundleSimulation:
    def __init__(self, T, dt, nerve, stimulusCurrent):
        self.T = T
        self.dt = dt
        self.timeLine = np.arange(0, mag((T+dt), ms), mag(dt, ms))
        self.nerve = nerve
        self.stimulusCurrent = stimulusCurrent

        # Create and place the axons
        print "Creating bundle of fibers..."
        for i in range(0, nerve["numFibers"]):
            succeeded = self.placeFiberInNerve()
            if not succeeded:
                break
        print "Placed " + str(len(nerve["fibers"])) + " fibers."

    def simulate(self, exciteCenterOnly=False):
        nerve = self.nerve
        stimulusCurrent = self.stimulusCurrent
        for t in range(1, len(self.timeLine)):
            if self.timeLine[t] % 1 == 0.0:
                print "Simulation Time: " + str(self.timeLine[t])

            for i, fiber in enumerate(nerve["fibers"]): # for each nerve fiber
                for k, axonNode in enumerate(fiber.axonNodes): # for each node in the fiber
                    distance = axonNode.distance
                    effectiveCurrent = getCurrent(t*self.dt, stimulusCurrent["magnitude"])

                    lastStep = len(axonNode.Vm) - 1

                    # initialize the left and right nodes to have the same values as the current node.
                    # this is so that they don't become current sinks (Boucher 2012 eq. 17 and 18)
                    leftNode  = {"V": axonNode.Vm[lastStep], "d": axonNode.distance, "n": axonNode.index}
                    rightNode = {"V": axonNode.Vm[lastStep], "d": axonNode.distance, "n": axonNode.index}

                    if (k-1) > -1:
                        node = fiber.axonNodes[k-1]
                        leftNode["V"] = node.Vm[lastStep]
                        leftNode["d"] = node.distance
                        leftNode["n"] = node.index

                    if (k+1) < len(fiber.axonNodes):
                        node = fiber.axonNodes[k+1]
                        rightNode["V"] = node.Vm[lastStep]
                        rightNode["d"] = node.distance
                        rightNode["n"] = node.index

                    # step the current axon forward IN TIIIME ♪♪
                    # print "Stepping axon #" + str(k) + " in fiber #" + str(i)
                    log.info("========")
                    log.infoVar(t*self.dt, 'time')
                    log.infoVar(axonNode.index, 'node')
                    log.infoVar(leftNode["V"], "leftVoltage")
                    log.infoVar(rightNode["V"], "rightVoltage")
                    axonNode.step(effectiveCurrent, leftNode, rightNode, self.dt, exciteCenterOnly)

    def dumpJSON(self, filename):
        import json
        outputFile = open(filename, 'w')

        data = {}
        data["timeLine"] = self.timeLine.tolist()
        data["stimulusCurrent"] = {
            "magnitude": self.stimulusCurrent["magnitude"],
            "x": self.stimulusCurrent["x"],
            "y": self.stimulusCurrent["y"],
            "z": self.stimulusCurrent["z"]
        }

        data["fibers"] = [] # an array of objects
        for j, fiber in enumerate(self.nerve["fibers"]):
            data["fibers"].append({
                "x": fiber.x,
                "y": fiber.y,
                "diameter": fiber.diameter,
                "internodalLength": fiber.internodalLength,
                "axonalLength": fiber.axonalLength,
                "nodes": []
            })

            for k, node in enumerate(fiber.axonNodes):
                Vm, m, h, n, mLS, hLS = node.Vm, node.m, node.h, node.n, node.mLS, node.hLS

                iNaSol = [node.sodiumCurrent(Vm[i], m[i], h[i], mLS[i], hLS[i]) for i in range(0, len(self.timeLine))]
                iKSol  = [node.potassiumCurrent(Vm[i], n[i]) for i in range(0, len(self.timeLine))]
                iLSol  = [node.leakageCurrent(Vm[i]) for i in range(0, len(self.timeLine))]
                curr = [getCurrent(t, self.stimulusCurrent["magnitude"]) for t in self.timeLine]
                extPotentialSol = [node.extV(curr[i], node.distance) for i, _ in enumerate(self.timeLine)]

                data["fibers"][j]["nodes"].append({
                    "index": node.index,
                    "distanceFromStimulus": node.distance,
                    "damagedChannels": node.damagedChannels,
                    "voltage": Vm,
                    "m": m,
                    "h": h,
                    "n": n,
                    "mLS": mLS,
                    "hLS": hLS,
                    "iNa": iNaSol,
                    "iK": iKSol,
                    "iL": iLSol,
                    "externalPotential": extPotentialSol
                })

        print "Writing file " + outputFile.name + "..."
        outputFile.write(json.dumps(data))
        print "Done."

    # Places a nerve fiber and makes sure it doesn't overlap with any other nerve fibers in the nerve bundle
    # returns true if it succeeds, false otherwise
    def placeFiberInNerve(self, maxAttempts = 1000):
        nerve = self.nerve
        current = self.stimulusCurrent
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

            return NerveFiber(x, y, diameter, nerve["numNodes"], nerve["axonalLength"], nerve["damageMap"], current)

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

# represents a square wave current strimulus
def getCurrent(t, current, tPulseStart=0*ms, pulseWidth=550*ms):
    if tPulseStart <= t <= (tPulseStart + pulseWidth):
        return current
    else:
        return 0.0*mA

def getDistance(x1, y1, z1, x2, y2, z2):
    x = abs(x1 - x2)
    y = abs(y1 - y2)
    z = abs(z1 - z2)
    return np.sqrt(mag((x**2 + y**2 + z**2), (cm*cm))) * (cm)

##########################
# Some Plotting Functions
##########################
# Checks if the destination directory exists and creates it if it doesn't.
# Use instead of `saveFigure`.
# I only tested it with paths/that/look/like/this.jpg
def saveFigure(path):
    folders = path.split('/')[:len(path.split('/'))-1]
    folderPath = "/".join(folders)
    if not os.path.isdir(folderPath):
        os.makedirs(folderPath)
    pylab.savefig(path)

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

    pylab.scatter(x, y, color='r', marker='o', s=10)
    pylab.ylabel('y (cm)')
    pylab.xlabel('x (cm)')
    pylab.axis('equal')
    for i, fiber in enumerate(nerve["fibers"]):
        circle = pylab.Circle((x[i],y[i]), mag(fiber.diameter/(2.0), cm), alpha=0.5)
        pylab.gca().add_artist(circle)

    if plotStimulusPos:
        pylab.axhline(y = mag(stimulusCurrent["y"], cm), color='k', linestyle='--')
        pylab.text(mag(stimulusCurrent["x"], cm), mag(stimulusCurrent["y"], cm), "stimulus")
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
    saveFigure("graphs/axons/axon" + str(node.index) + "fiber" + str(fiberNum) + ".jpg")
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

def plotInfoOfNodes(simulation, plotStimulus=True, exciteCenterOnly=False):
    curr = []
    for j in range(0, len(simulation.timeLine)):
        curr.append(getCurrent(j*simulation.dt, simulation.stimulusCurrent["magnitude"]))

    for fiberNum, fiber in enumerate(simulation.nerve["fibers"]):
        for k, node in enumerate(fiber.axonNodes):
            print "plotting axon #" + str(node.index) + " in fiber #" + str(fiberNum) + "..."
            node.plotCurrentsVoltagesAndGates(simulation.timeLine, curr, fiberNum, plotStimulus, exciteCenterOnly)

def plotCompoundPotential(simulation, n=0):
    compoundPotential = [0 for i in range(0, int(simulation.T/simulation.dt) + 1)]

    arrayIndex = n + int(simulation.nerve["numNodes"]/2)
    # add up the nodes in each fiber with the same n. These nodes all have the same index,
    # but are not necessarily next to each other. They should be close though
    for fiber in simulation.nerve['fibers']:
        node = fiber.axonNodes[arrayIndex]
        for k, v in enumerate(node.Vm):
            compoundPotential[k] += v

    pylab.figure()
    pylab.plot(simulation.timeLine, compoundPotential)
    pylab.xlabel('Time (ms)')
    pylab.ylabel('Voltage (mV)')
    pylab.title('Sum of all Action Potentials at Node n = ' + str(n))
    pylab.grid()
    saveFigure("graphs/sumOfPotentials.jpg")
    pylab.close()

##############
# Start Script
##############
def main():
    log.logLevel = log.ERROR
    # log.logLevel = log.INFO

    # parse arguments
    parser = argparse.ArgumentParser(description="Nerve Fiber Loss Simulation Tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--numfibers', default=1, type=int, help='the number of neuron fibers to simulate')
    parser.add_argument('--numnodes', default=11, type=int, help='the number of axon nodes in each fiber to simulate')
    parser.add_argument('--time', default=55, type=float, help='the amount of time to simulate in ms')
    parser.add_argument('--dt', default=0.025, type=float, help='the time step to use')

    args = parser.parse_args()

    # Current Stimulus
    # threshold is ~12uA/cm^2 or 2.9e-6uA
    stimulusCurrent = {
        "magnitude" : -2.59 *uA,    # uA. the current applied at the surface
        "x"         : 0   *cm,    # cm
        "y"         : 0.3 *cm,    # cm
        "z"         : 0   *cm     # cm
    }

    # the nerve is a bundle of nerve fibers. Nerve fibers are rods of connected axons.
    nerve = {
        "numFibers"    : args.numfibers,
        "numNodes"     : args.numnodes,    # the number of axon nodes each fiber has. Should be an odd number.
        "fibers"       : [],
        "radius"       : 0.0    *cm, # cm
        "x"            : 0.0    *cm, # cm
        "y"            : 0.0    *cm, # cm
        "z"            : 0.0    *cm, # cm
        "minFiberDiam" : 0.0019 *cm, # cm
        "maxFiberDiam" : 0.0021 *cm, # cm
        "axonalLength" : 2.5e-4 *cm, # cm
        "damageMap"       : {           # a map of node indices to the proportion of damaged channels on that node.
            11 : 1.0,
        }
    }

    # plotNodePositions()
    # plotCrossSectionPositions()

    T    = args.time*ms    # ms
    dt   = args.dt*ms # ms
    centerOnly = False

    simulation = NerveBundleSimulation(T, dt, nerve, stimulusCurrent)
    print "Starting simulation... Length =", str(T), "ms. Step =", str(dt), "ms.", "Points:", str(T/dt)
    print "numFibers =", str(nerve["numFibers"]), "numNodes =", str(nerve["numNodes"])
    simulation.simulate(exciteCenterOnly=centerOnly) # modifies `nerve`
    plotInfoOfNodes(simulation, plotStimulus=True, exciteCenterOnly=centerOnly)
    plotCompoundPotential(simulation, n=11)
    # simulation.dumpJSON('data.json')

    print "Done."

if  __name__ =='__main__':
    main()
