import math
from fractions import Fraction
import pandas as pd
from primes import primes
import itertools
from dataclasses import dataclass
from typing import Callable

primes_max = 10000

As = [
    # [0],
    # [0,2],
    # [0,6],
    # [0,10],
    # [0,14],
    # [0,30],
    # [0,2,6],
    # [0,4,10],
    # [0,4,28],
    # [0,2,6,8],
    # [0,18,30,60],
    # [0,60,120,180],
    # [0,2,6,8,30],

    [0,2],        
    [0,10],       
    [0,6],        
    [0,2,6,8],    
    [0,8,12,30],  
    [0,6,12,18],  
    [0,2,6,8,30], 
    [0,18,30,60], 
    [0,60,120,18],
    ]

# Compute residue counts
def r_gen(A): # u for primes
    def out(p):
        A_mod_p = set()
        for k in A:
            A_mod_p.add( k % p )
        return len(A_mod_p)
    return out

# Find modulated primes for each uA
def find_mods(A, rA):
    mods = []
    for p in primes[:primes_max]:
        if rA(p) < len(A):
            mods.append(p)
    return mods

# Compute M_A(1)
def compute_M_A(A, r, mods):
    c = len(A)
    return math.prod( (1-Fraction(1,p)) / (1-Fraction(r(p),p*c)) for p in mods )

# Compute N_A(1)
def compute_N_A(A, r):
    c = len(A)
    return math.prod( (1-1/p)**c / (1-r(p)/p) for p in primes[:primes_max] )

@dataclass
class ktuple:
    A: list
    r: Callable[[int], int]
    mods: list
    M: float
    N: float

kts = []
for A in As:
    r = r_gen(A)
    if any( r(p) == p for p in primes[:100]):
        continue
    mods = find_mods(A, r)
    M = compute_M_A(A, r, mods)
    N = compute_N_A(A, r)
    kts.append( ktuple(A, r, mods, M, N) )

for k in kts:
    print( round(float(k.M),8), round(k.N, 8), round(1.0/k.N, 8), k.A, k.mods )

# Print everything 
# pd.set_option("display.float_format", "{:.8f}".format)
# print()
# print(pd.DataFrame({
#     "k-tuple"   : kts,
#     "modulates" : mods,
#     "M_A frac"  : M_As,
#     "M_A dec"   : [float(M) for M in M_As],
#     "N_A"       : [float(N) for N in N_As],
# }))
# print()  