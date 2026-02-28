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
    [0,60,120,180],
    ]

# Compute residue counts
def r_gen(A): # u for primes
    def out(p):
        A_mod_p = set()
        for k in A:
            A_mod_p.add( k % p )
        return len(A_mod_p)
    return out

# Takes pairs of residues and moduli, finds a simultaneous solution 
def CRT( pairs ):
    def mod_inverse(a, m):
        # Since m is prime, we can use Fermat's little theorem: a^(m-1) ≡ 1 (mod m)
        # So a^(-1) ≡ a^(m-2) (mod m)
        return pow(a, m - 2, m)
    M = math.prod(mi for _,mi in pairs)
    result = 0
    for ai, mi in pairs:
        Mi = M // mi
        yi = mod_inverse(Mi, mi)
        result += ai * Mi * yi
    return result % M

# Theorem 1.14 - algorithm for finding k-tuple with given residues
def find_r_rep( r, N ):
    # Get M
    M = 0
    for p in primes[:N]:
        if M < r(p): M = r(p)
    A_i = [0]
    for i in range(2, M+1):
        # Now, find r for current Ai and check at which primes it's smaller than our desired r
        ri = r_gen( A_i )
        Q_i = [p for p in primes[:N] if ri(p) < r(p)]
        
        # CRT gives a k with 0 <= k < pN# that solves 
        #   k congruent i (mod q) for all q in Q_i
        #   and
        #   k congruent 0 (mod p_i) for p_i in (Primes - Q_i) cap [0,pN]
        CRT_package = []
        for p in primes[:N]:
            if any( p == q for q in Q_i ):
                CRT_package.append([i,p])
            else:
                CRT_package.append([0,p])                
        k = CRT( CRT_package )

        A_i.append( k )
    A_i.sort()
    return A_i

# Find modulated primes for each uA
def find_mods(A_size, rA, A_max_diff):
    mods = []
    for p in primes[:primes_max]:
        if p > A_max_diff: break
        if rA(p) < A_size:
            mods.append(p)
    return mods

# Compute M_A(1)
def compute_M_A(A_size, r, mods):
    k = A_size
    return math.prod( (1-Fraction(1,p)) / (1-Fraction(r(p),p*k)) for p in mods )

# Compute N_A(1)
def compute_N_A(A_size, r):
    k = A_size
    return math.prod( (1-1/p)**k / (1-r(p)/p) for p in primes[:primes_max] )

@dataclass
class ktuple:
    # A: list
    A_len: int
    A_rep: list
    r: Callable[[int], int]
    mods: list
    M: float
    N: float

# kts = []
# for A in As:
#     r = r_gen(A)
#     if any( r(p) == p for p in primes[:100]):
#         continue
#     mods = find_mods(len(A), r)
#     M = compute_M_A(len(A), r, mods)
#     N = compute_N_A(len(A), r)
#     kts.append( ktuple(A, r, mods, M, N) )

# Brute force search for generating tuple
# Currently only works for k=3
def search_for_generating_tuple( mods, residues, bound=1000 ):
    for a1 in range(2,bound):
        for a2 in range(1,a1):
            A = [0,a2*2,a1*2]
            r = r_gen(A)

            # First, make sure all the mods match all the residues
            if any( r(mod) != residue for mod, residue in zip( mods, residues ) ):
                continue
            
            # Now make sure there are no new mods
            new_mods = find_mods( 3, r, A[2] )
            if not mods == new_mods:
                continue

            # If both checks are passed, we found a generating tuple
            return A
    return None

bad_units = []
print()
# Check different combos of modded primes
good = []
bad = []
for k3 in range(2,40+1):
    # for k2 in range(1,k3):
        # for k1 in range(0,k2):
            ks = [0,1,k3] # 2 and 3 must be modulated for admissibility
            ps = [primes[k] for k in ks]
            # And then for those combos try different residue counts
            for r3 in range(1,min(ps[2],3)):
                for r2 in range(1,min(ps[1],3)):
                    for r1 in range(1,min(ps[0],3)):
                        rs = [r1,r2,r3]
                        # Now we have one unit function!
                        # Lets try to find a k-tuple that generates it
                        kt = search_for_generating_tuple( ps, rs )
                        # def r(p):
                        #     if p == ps[0]: return rs[0]
                        #     elif p == ps[1]: return rs[1]
                        #     elif p == ps[2]: return rs[2]
                        #     else: return 3
                        if kt == None: 
                            bad.append( [ps,rs] )
                        else: 
                            good.append( [ps,rs] ) 
print()

def disp( l ):
    # Ms = [float( compute_M_A(3, r, ps) ) for ps,rs in l]
    data = []
    for ps,rs in l:
        def r(p):
            if p == ps[0]: return rs[0]
            elif p == ps[1]: return rs[1]
            elif p == ps[2]: return rs[2]
            else: return 3
        data.append([compute_M_A(3, r, ps), compute_N_A(3, r),ps,rs])

    # Ms = list( set(Ms) )
    data.sort(key=lambda d: d[0])
    for datum in data: print( datum[0], datum[1], datum[2][2], datum[3] )

print( "BAD:" )
disp( bad )
# print( "GOOD:" )
# disp( good )

# A_len = 3
# def r(p):
#     if p == 2: return 1
#     elif p == 11: return 2
#     else: return 3
# # mods = find_mods( A_len, r )
# rep = find_r_rep( r, A_len )
# print( [r(2),r(3),r(5),r(7),r(11),r(13)], rep )


# kts.append( ktuple(A_len, rep, r, mods,
#     compute_M_A(A_len, r, mods),
#     compute_N_A(A_len, r) ) )

# print()
# for k in kts[:1]:
#     print( round(float(k.M),8), round(k.N, 8), round(1.0/k.N, 8), k.A_rep, k.mods )
# print()

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