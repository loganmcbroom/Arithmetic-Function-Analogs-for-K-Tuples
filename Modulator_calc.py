import math
from fractions import Fraction
import pandas as pd
from primes import primes
import itertools
from dataclasses import dataclass
from typing import Callable

primes_max = 1000000

# We need to create a bunch of kts which are admissible and non-residually equivalent
# We can do 4 per kt and make them all multiples of 6 to guarantee adm.

# kts = [
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
    # ]

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
for A in itertools.combinations(range(1,6),3):
    A = list(A)
    A.insert(0,0)
    A = [k*6 for k in A]

    def gen_ktuple( A ):
        r = r_gen(A)
        if any( r(p) == p for p in primes[:100]):
            return None
        mods = find_mods(A, r)
        M = compute_M_A(A, r, mods)
        N = round( compute_N_A(A, r), 9 )
        return ktuple(A, r, mods, M, N)

    kts.append( gen_ktuple( A.copy() ) ) 
    A[1] += 2
    kts.append( gen_ktuple( A ) ) 

def dedupe_by_key(seq, key):
    seen = set()
    result = []
    for item in seq:
        k = key(item)
        if k not in seen:
            seen.add(k)
            result.append(item)
    return result
from operator import attrgetter
kts = dedupe_by_key(kts, key=attrgetter("N"))

kts.sort(key = lambda x: x.M)

def is_sorted(seq, key=lambda x: x):
    return all(key(a) <= key(b) for a, b in zip(seq, seq[1:]))

for k in kts:
    print( float(k.M), k.N, k.A, k.mods )

print( is_sorted( kts, key=lambda x: x.N ) )

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