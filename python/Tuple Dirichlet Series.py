from base import *
from math import exp, log, prod
# import numpy as np
# import matplotlib.pyplot as plt

def TDS( A, f ):
    def out(s):
        total = 0
        for n in range(1,1000):
            total += prod(f(a+n) for a in A) / exp(s*prod(log(a+n) for a in A))
        return total
    return out
    # return np.vectorize( out )

'''
Searching for TDS product collapsing in some way, but it seems to never collapse?
As in, the TDS product has no terms in the same convolution bin.
This seems really unlikely, maybe FE-type conjecture at play?
'''
for tb in range(5,50):
    ta = 2
    def f(a,b):
        return a*(a+2) * b*(b+2)
    test = f(ta,tb)
    for a in range(2,500):
        for b in range(a,500):
            if test == f(a,b):
                if a != ta or b != tb:
                    print( a,b, ta,tb )
# A = [0,2,6]
# s = 2
# La = conv(invert([1]*bound), eval(log)) # Von Mangoldt
# for i in range(bound): 
#     if abs(La[i]) < 1e-14: La[i] = 0.0
# La_A_series = TDS(A, lambda n: La[n])(s)
# zp_over_z = TDS(A, lambda n: log(n) )(s) / TDS(A, lambda _: 1)(s)
# print( La[:20] )