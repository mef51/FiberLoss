#!/usr/bin/python

from unum.units import * # this import is slow

c = 0*m
its = 100000
for i in range(0, its):
    c += 1*m
