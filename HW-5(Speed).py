import math
import random
import statistics
import numpy
from timeit import default_timer as timer
from numba import jit


def v_asian_sample(T, r, K, s, X0, m):
    xhat = 0.0
    X = X0
    D = T / m
    for i in range(m):
        X *= math.exp(random.normalvariate((r-s**2/2)*D, s*D**0.5))
        xhat += X
    return math.exp(-r*T)*max(xhat/m - K, 0)

@jit
def v_asian_sample_jit(T, r, K, s, X0, m):
    xhat = 0.0
    X = X0
    D = T / m
    for i in range(m):
        X *= math.exp(random.normalvariate((r-s**2/2)*D, s*D**0.5))
        xhat += X
    return math.exp(-r*T)*max(xhat/m - K, 0)

def v_asian_sample_vec(T, r, K, s, X0, m):
    D = T / m
    X = numpy.random.normal((r-s**2/2)*D, s*D**0.5, m)
    return math.exp(-r*T)*max(numpy.mean(numpy.exp(numpy.cumsum(X)))*X0 - K, 0)


def v_asian(T, r, K, s, Xo, m, n, fun):
    return statistics.mean([fun(T, r, K, s, Xo, m) for i in range(n)])

for f in [v_asian_sample, v_asian_sample_jit, v_asian_sample_vec]:
    for i in range(2):
        start = timer()
        answer = v_asian(1.0, 0.05, 55.0, 0.3, 50, 100, 5000, f)
        print(f, answer)
        duration = timer() - start
        print(f, "\n  ", duration, " seconds")