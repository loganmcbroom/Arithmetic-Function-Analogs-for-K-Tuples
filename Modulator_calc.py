import math
from fractions import Fraction
import pandas as pd
from primes import primes

primes_max = 1000000

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
    [0,18,30,60],
    [0,60,120,180],
    [0,2,6,8,30],
    
    # [0,2,6,8],
    # [0,6,12,18],
    # [0,2,6,30],
    # [0,120,122,126],
    # [0,8,120,180],
    # [0,12,22,42],
    # [0,6,30,42],

    ]

# Compute residue counts
def r_gen(A): # u for primes
    def out(p):
        A_mod_p = set()
        for k in A:
            A_mod_p.add( k % p )
        return len(A_mod_p)
    return out
rAs = [r_gen(A) for A in kts]

# Find modulated primes for each uA
def find_mods(A, rA):
    mods = []
    for p in primes[:primes_max]:
        if rA(p) < len(A):
            mods.append(p)
    return mods
mods = [find_mods(kts[i], rAs[i]) for i in range(len(kts))]

# Compute M_A(1)
def compute_M_A(A, rA, mod):
    c = len(A)
    return math.prod( (1-Fraction(1,p)) / (1-Fraction(rA(p),p*c)) for p in mod )
M_As = [compute_M_A(kts[i], rAs[i], mods[i]) for i in range(len(kts))]

# Compute N_A(1)
def compute_N_A(A, rA):
    c = len(A)
    return math.prod( (1-1/p)**c / (1-rA(p)/p) for p in primes[:primes_max] )
N_As = [compute_N_A(kts[i], rAs[i]) for i in range(len(kts))]


# Print everything 
pd.set_option("display.float_format", "{:.8f}".format)
print()
print(pd.DataFrame({
    "k-tuple"   : kts,
    "modulates" : mods,
    "M_A frac"  : M_As,
    "M_A dec"   : [float(M) for M in M_As],
    "N_A"       : [float(N) for N in N_As],
}))
print()  