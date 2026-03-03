from base import *

"""
This file is for extensions that demand types of incidence
other than the usual complete containment.

We currently do everything by count,
but we can be more exact by demanding certain portions of the k-tuple
meet a condition, and the rest do not, e.g., for phi we can
ask for the 0 of [0,2] to be relatively prime and the 2 to not.

It seems that demanding less than full containment causes the 
associated u function to be eventually 0 on primes.
It seems like these u functions are still determined
in some way by residual collisions, but it isn't clear
to me what they are.
"""

# Generate phi extensions that meet a set incidence count demand,
# either exactly, or at least.
def ophi_A_gen(A, demand = 4, exact = False):
    if demand == -1: demand = len(A)
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

# This gives all the u extensions on primes
A = [0,2,6,12]
print(A)
print(2, 3, 5, 7, 11)  
for demand in range(1,len(A)+1):
    ophi_2 = ophi_A_gen(A, demand)
    print(2-ophi_2[2], 3-ophi_2[3], 5-ophi_2[5], 7-ophi_2[7], 11-ophi_2[11])   