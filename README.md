Fiber Loss in Nerves
====================

People who suffer from diseases that cause nerve fiber loss often don't realize what's happening because it happens so slowly.
Tests and experiments can be done where you measure how far a signal travels down a limb from a certain stimulus. The more damage there is, the less the signal travels.

For example, an electrode is placed near your elbow, and a voltmeter is placed further down the arm, around your wrist. The electrode stimulates some nerves in your elbow and the signal should propagate down your arm to the voltmeter, where we can measure a signal. If the nerves in your arm are damaged in some way, this signal will be weaker.

This simulation is inspired by a mathematical model that gives some clues as to the nature of this damage. The model seems to be able to explain many kinds of damage as the result of a single malfunction, namely, leaky ion channels in the cellular membrane of the nerve.

The goal of the simulation is to recreate results from the experiment described above using the mathematical model.

Status: Creates virtual nerve fibers, generates action potentials, and propagats to connected nodes
====================

* Creates a bundle of fibers, where each fiber is a group of connected nodes.
* Places a stimulus some distance away from the bundle.
* Can simulate the response of each node to the axon. Axons affect each other, and can potentially trigger propagating action potentials.
* Can plot various things about the system, including positions, membrane voltages, and forward and backward rates.

Papers
======
This work is mostly based off of the theory in these two papers:

* [Analysis of a Model for Excitation of Myelinated Nerve (Mcneal 1976)](papers/McNealModel.pdf)
* [Coupled left-shift of Nav channels (Boucher, Joos, Morris, 2012)](papers/BoucherCLS.pdf)

For more details on the science this is based off of or some related papers take a look at the [papers/](papers/) folder.
