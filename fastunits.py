# This file is meant to replace the import for the unum.units package.
# It redefines the units to be '1.0' so that simulations run faster.
# Using unum.units is great because it helps you catch errors and lets you experiment safely while developping.
# However this unit checking takes a long time and makes simulation runs much slower.
# Whenever a real simulation is to be run, just import this file instead of the unum.units package.
# This file is meant to keep your existing code running as it is (with all the units and stuff) and to speed up simulations
# when you're confident the units are right.

# Redefine Unum to do nothing but have the same interface
class Unum(object):
    def __init__(self, value):
        self.value = value

    def unit(cls, s, n):
        return 1.0
    unit = classmethod(unit)

V    = 1.0
mV   = 1.0
uV   = 1.0

F    = 1.0
mF   = 1.0
uF   = 1.0

m    = 1.0
mm   = 1.0
cm   = 1.0

ohm  = 1.0
mohm = 1.0

S    = 1.0
mS   = 1.0

A    = 1.0
mA   = 1.0
uA   = 1.0

s    = 1.0
ms   = 1.0


