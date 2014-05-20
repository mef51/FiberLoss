import numpy as np
import pylab


def plotL():
    Ds = np.arange(0, 3e-5, 0.1e-5)
    Ls = []
    for D in Ds:
        Ls.append(internodalLength(D))

    pylab.plot(Ds, Ls)
    pylab.title("Fibre Diameter vs. Internodal Length")
    pylab.show()

def integrateExample():
    f = lambda x: x**2
    dx = 0.00000001
    interval = np.arange(0, 2, dx)
    integralOfX = 0
    for i, xi in enumerate(interval):
        integralOfX += f(xi) * dx
    print integralOfX
