from base import *
from itertools import combinations 

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

Slept on that, you count how many of the residue classes
obtained a certain number of fills.
Does exact mode require exactly the number of fills?
"""

mu = invert([1]*bound)
N = range(bound)
phi = conv(N,mu)
def mu_gen(A):
    return pointwise(u_gen(A), mu)
def phi_gen(A):
    return conv(N, mu_gen(A))

# Generate phi extensions that meet a set incidence count demand,
# either exactly, or at least.
def ophi_A_gen(A, quota, exact):
    def out(n):
        total = 0
        for k in range(n):
            corrects = set()
            for a in A:
                if math.gcd(k+a,n) == 1:
                    corrects.add((k+a)%n)
            met = len(corrects)
            if exact:
                if quota == met: total += 1
            else:
                if quota <= met: total += 1
        return total
    return eval( out )

def composed(n,B): # Test if n is composed of B
    fs = factorization(n)
    return set(fs).issubset( B )

def ophi_A_smooth(A, B):
    def out(n):
        total = 0
        for k in range(n):
            if all( composed(math.gcd(k+a,n),B) for a in A):
                total += 1
        return total
    return eval( out )

# Generate an exact quota phi in terms of full quotas
def synthesize_exact_quota(A,q):
    ps = [ [phi_gen(s) for s in combinations(A,r)] for r in range(1,len(A)+1) ]
    return [int( sum( 
        (-1)**(k-q) * math.comb(k,q) * sum(p[i] for p in ps[k-1] ) 
        for k in range(q, len(A)+1))
    ) for i in range(bound)]
    
# Generate a minimum quota phi in terms of full quotas
def synthesize_partial_quota(A,q):
    ps = [ [phi_gen(s) for s in combinations(A,r)] for r in range(1,len(A)+1) ]
    return [int( sum( 
            (-1)**(k-q) * math.comb(max(k-1,0),k-q) * sum(p[i] for p in ps[k-1] ) 
            for k in range(q, len(A)+1)) 
        ) for i in range(bound)] 


A = [0,2,6,12,30,60]

ophi = ophi_A_smooth( A, set([2,3]) )
u_smooth = conv(invert(ophi),N)
for i in range(bound): u_smooth[i] = int(u_smooth[i])
print(ophi[:20])
print(u_smooth[:20])

# Tests for forcing quotad phis to be multiplicative
# ophi = ophi_A_gen(A, q)
# print(ophi[:30])

# def r(p):
#     return p - ophi[p]

# def gen_forced_mult_phi(r):
#     def out(n):
#         ps = [p for p in primes[:100] if n%p == 0]
#         return int( n * math.prod( 1.0 - r(p)/p for p in ps ) )
#     return eval(out)

# p = gen_forced_mult_phi(r)
# print(p[:30])


