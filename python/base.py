import math
from primes import primes
from fractions import Fraction

bound = 2**10+1 # Computation bound
primes_max = 10000 # Largest prime list index to use

def factorization(n):
    factors = {}
    for p in primes:
        if p > n: break
        while n % p == 0:
            if p in factors: factors[p] += 1
            else: factors.update({p:1})
            n = n // p
        if n == 1: break
    return factors

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

def eval( f ):
    return [f(i) if i > 0 else 0 for i in range(0,bound)]

Id = [0]*(bound-1)
Id.insert( 1, 1 )

def conv_power( f, n, g = Id ):
    if n == 0: return f
    else: return conv_power( conv( f, g ), n-1, f )

# u = [1]*bound
# mu = invert( u )
# phi = conv( mu, range(bound) )
# La = conv(mu, eval(lambda n: math.log(n))) # Von Mangoldt
# for i in range(bound):
#     if abs(La[i]) < 1e-14: La[i] = 0.0

# Liouville 
# def is_square(n):
#     return int(math.sqrt(n)) ** 2 == n
# sq = [1 if is_square(i) else 0 for i in range(bound)]
# la = conv(mu,sq) 

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
        return math.prod(r(p)**e for [p,e] in factor(n).items())
    return eval(out)

def r(p, A): # u for primes
    A_mod_p = set()
    for k in A:
        A_mod_p.add( k % p )
    return len(A_mod_p)

# Compute residue counts
def r_gen(A): # u for primes
    def out(p):
        A_mod_p = set()
        for k in A:
            A_mod_p.add( k % p )
        return len(A_mod_p)
    return out

def u_gen(A):
    return unit_gen( A, r_gen(A) )

def w_gen(A):
    return unit_gen( A, lambda p : r_gen(A)(p)/len(A) )

def mu_gen(A):
    return pointwise(u_gen(A), mu)

def phi_gen(A):
    return conv(N, mu_gen(A))

# Find modulated primes for each uA
def find_mods(A):
    r = r_gen(A)
    A_max_diff = A[-1] - A[0]
    mods = []
    for p in primes[:primes_max]:
        if p > A_max_diff: break
        if r(p) < len(A):
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

# Compute the singular series for A
def singular_series_of_A(A):
    return 1.0/compute_N_A(len(A), r_gen(A))
