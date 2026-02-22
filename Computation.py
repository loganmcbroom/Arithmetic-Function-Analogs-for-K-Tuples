import math
import cmath

bound = 3**8+1

def divisors(n):
    n = abs(n)
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = 1
    divs = [1]
    for p, e in f.items():
        current = []
        for d in divs:
            for k in range(1, e + 1):
                current.append(d * p**k)
        divs += current
    return sorted(divs)
def prod_bp(n):
    return [(d,n//d) for d in divisors(n)]

# Dirichlet product
def conv( f, g ):
    return [ sum( f[a]*g[b] for a,b in prod_bp(n) ) for n in range(bound) ]

# Pointwise product
def pointwise( f, g ):
    return [ f[n]*g[n] for n in range(bound) ]

def U( f ):
    return [ sum( f[a] for a,b in prod_bp(n) if b == n ) for n in range(bound) ]

# Dirichlet inverse
def invert( f ):
    f_inv = [None]*bound
    U_f = U( f )
    if not all( u != 0 for u in U_f ):
        print( "Could not invert function" )
        return f_inv
    f_inv[0] = 0
    f_inv[1] = 1/f[1]
    for n in range(2, bound):
        f_inv[n] = -sum([f[a]*f_inv[b] for a,b in prod_bp(n) if b < n ]) / U_f[n]
    return f_inv

# Base functions ========================================================
def eval( f ):
    return [f(i) if i > 0 else 0 for i in range(0,bound)]

u = [1]*bound
mu = invert( u )
N = range(bound)
phi = conv( mu, N )

from collections import defaultdict
def unit_gen(A, r):
    def factor(n: int) -> dict[int, int]:
        factors = defaultdict(int)
        while n % 2 == 0:
            factors[2] += 1
            n //= 2
        p = 3
        while p * p <= n:
            while n % p == 0:
                factors[p] += 1
                n //= p
            p += 2
        if n > 1:
            factors[n] += 1
        return dict(factors)
    def out(n):
        if n < 1: return 0
        return math.prod(r(p, A)**e for [p,e] in factor(n).items())
    return eval(out)
def r(p, A): # u for primes
    A_mod_p = set()
    for k in A:
        A_mod_p.add( k % p )
    return len(A_mod_p)
def u_gen(A):
    return unit_gen( A, r )
def w_gen(A):
    return unit_gen( A, lambda p : r(p)/len(A) )

def mu_gen(A):
    return pointwise(u_gen(A), mu)

def phi_gen(A):
    return conv(N, mu_gen(A))

# A = set([0,2,6])
# u_A = u_gen(A)
# w_A = w_gen(A)
# mu_A = mu_gen(A)
# phi_A = phi_gen(A)
# def is_square(n):
#     return int(math.sqrt(n)) ** 2 == n
# sq = [1 if is_square(i) else 0 for i in range(bound)]
# la = conv(sq,mu)

# Dirichlet character extension testing
# w = cmath.exp((2*math.pi / 6.0)*1j)
# w2 = w**2
# chis = [
#     [0,1,0,  1,0,  1,0,0,0,  1,0,  1,0,  1]*64,
#     [0,1,0,  w,0,-w2,0,0,0, w2,0, -w,0, -1]*64,
#     [0,1,0, w2,0, -w,0,0,0, -w,0, w2,0,  1]*64,
#     [0,1,0, -1,0, -1,0,0,0,  1,0,  1,0, -1]*64,
#     [0,1,0, -w,0, w2,0,0,0, w2,0, -w,0,  1]*64,
#     [0,1,0,-w2,0,  w,0,0,0, -w,0, w2,0, -1]*64,
#     ]
# chi_As = [pointwise(chi,u_A) for chi in chis]

# Omega testing
# def Omega(n):
#     O = 0
#     d = 2
#     while d <= n:
#         while n % d == 0:
#             O += 1
#             n //= d
#         d += 1
#     return O
# c_Omega = eval( lambda n : len(A)**Omega(n) )

# f = eval(lambda n : c_Omega[n]/u_A[n])
# for n in range(1,20):
#     print( c_Omega[n], end = ' ' )

# Experimental ======================================================

# lambda_A, based on L(s,la) = Z(2s)/Z(s)
# Instead of Z(2s), though, we use L(2s,u_A^2)
# That also gives L(s,la_A) = Prod 1/(1+u_A(p)p^(-s))
# la_A = conv(pointwise( sq, u_A ), mu_A)
# I proved this is equal to the following:
# la_A = pointwise( la, u_A )

# print( conv(w_A, mu) )


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
    ]

muAs = [mu_gen(A) for A in kts]

log = eval(lambda n : math.log(n))
LAs = [conv(mu_A, log) for mu_A in muAs]
from fractions import Fraction
for LA in LAs:
    for i in range(len(LA)):
        LA[i] = round( math.exp(LA[i]), 3 )

from primes import primes
def r_gen(A): # u for primes
    def out(p):
        A_mod_p = set()
        for k in A:
            A_mod_p.add( k % p )
        return len(A_mod_p)
    return out
rAs = [r_gen(A) for A in kts]
def find_mods(A, rA):
    mods = []
    for p in primes[:100]:
        if rA(p) < len(A):
            mods.append(p)
    return mods
mods = [find_mods(kts[i], rAs[i]) for i in range(len(kts))]

import pandas as pd
pd.set_option("display.float_format", "{:.3f}".format)
pd.set_option("display.max_colwidth", 200 )
print()
print(pd.DataFrame(LAs
    # "k-tuple"   : kts,
    # "modulates" : mods,
    # "Lambda_A"  : LAs
).T)
print() 