from base import *

"""
This file is for extensions that demand different types of incidence, rather
than the usual complete containment.

We currently do everything by count,
but we can be more exact by demanding certain portions of the k-tuple
meet a condition, and the rest do not, e.g., for phi we can
ask for the 0 of [0,2] to be relatively prime and the 2 to not.
"""

def ophi_A_gen(A, demand = 1, exact = False):
    def out(n):
        total = 0
        for k in range(n):
            met = 0
            for a in A:
                if math.gcd(k+a,n) == 1:
                    met += 1
            if exact:
                if demand == met:
                    total += 1
            else:
                if demand <= met:
                    total += 1
        return total
    return eval( out )

ophi_2 = ophi_A_gen([0,6])
print(ophi_2[:30])   