import Computation
import math
from fractions import Fraction
import pandas as pd

primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]

kts = [
    [0],
    [0,2],
    [0,6],
    [0,10],
    [0,14],
    [0,30],
    [0,2,6],
    [0,4,10],
    [0,4,28],
    [0,2,6,8],
    [0,2,6,30],
    ]

# Compute u_As, we only need primes but that's alright
uAs = [Computation.u_gen(A) for A in kts]

# Find modulated primes for each uA
def find_mods(A, uA):
    mods = []
    for p in primes:
        if uA[p] < len(A):
            mods.append(p)
    return mods
mods = [find_mods(kts[i], uAs[i]) for i in range(len(kts))]

# Compute M_A(1)
def compute_M_A(A, uA, mod):
    c = len(A)
    return math.prod( (1-Fraction(1,p)) / (1-Fraction(uA[p],p*c)) for p in mod )
M_As = [compute_M_A(kts[i], uAs[i], mods[i]) for i in range(len(kts))]

# Compute N_A(1)
def compute_N_A(A, uA, mod):
    c = len(A)
    return math.prod( (1-Fraction(1,p))**c / (1-Fraction(uA[p],p)) for p in mod )
N_As = [compute_N_A(kts[i], uAs[i], mods[i]) for i in range(len(kts))]


# Print everything 
print()
print(pd.DataFrame({
    "k-tuple"   : kts,
    "modulates" : mods,
    "M_A"       : M_As,
    "M_A dec"   : [float(M) for M in M_As],
    "N_A"       : N_As,
    "N_A dec"   : [float(N) for N in N_As],
}))
print()