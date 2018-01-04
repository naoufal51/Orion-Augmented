from numpy import *
def nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n
