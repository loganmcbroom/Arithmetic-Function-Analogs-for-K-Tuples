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
phi = to_int( conv(N,mu) )
def mu_gen(A):
    return pointwise(u_gen(A), mu)
def phi_gen(A):
    return to_int( conv(N, mu_gen(A)) )

def upsilon_gen(A, quota, exact):
    def out(n):
        total = 0
        for k in range(n):
            met = 0
            for a in A:
                if math.gcd(k+a,n) == 1:
                    met += 1
            if exact:
                if quota == met: total += 1
            else:
                if quota <= met: total += 1
        return total
    return eval( out )

def rho_gen(A, quota, exact):
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

# Generate an exact quota totient in terms of full quotas
def synthesize_ups_exact(A,q):
    ps = [ [phi_gen(s) for s in combinations(A,r)] for r in range(1,len(A)+1) ]
    return [int( sum( 
        (-1)**(k-q) * math.comb(k,q) * sum(p[i] for p in ps[k-1] ) 
        for k in range(q, len(A)+1))
    ) for i in range(bound)]
    
# Generate a minimum quota totient in terms of full quotas
def synthesize_ups_minimum(A,q):
    ps = [ [phi_gen(s) for s in combinations(A,r)] for r in range(1,len(A)+1) ]
    return [int( sum( 
            (-1)**(k-q) * math.comb(max(k-1,0),k-q) * sum(p[i] for p in ps[k-1] ) 
            for k in range(q, len(A)+1)) 
        ) for i in range(bound)] 

def synthesize_rho_exact(A,q):
    out = []
    for i in range(bound):
        Ai = wrap(A,i)
        ps = [ [phi_gen(s) for s in combinations(Ai,r)] for r in range(1,len(Ai)+1) ]
        val = 0
        for k in range(q, len(Ai)+1): 
            val += (-1)**(k-q) * math.comb(k,q) * sum(p[i] for p in ps[k-1] ) 
        out.append(int( val ) )
    return out

def synthesize_rho_minimum(A,q):
    out = []
    for i in range(bound):
        Ai = wrap(A,i)
        ps = [ [phi_gen(s) for s in combinations(Ai,r)] for r in range(1,len(Ai)+1) ]
        val = 0
        for k in range(q, len(Ai)+1): 
            val += (-1)**(k-q) * math.comb(max(k-1,0),k-q) * sum(p[i] for p in ps[k-1] ) 
        out.append(int( val ) )
    return out

# ===================================================================================================
# ===================================================================================================
# ===================================================================================================


# Testing for jump functions between totient levels
# A = [0,2,6]
# phi_A = phi_gen(A)

# p0 = upsilon_gen(A,0,False)
# p1 = upsilon_gen(A,1,False)
# p2 = upsilon_gen(A,2,False)
# p3 = upsilon_gen(A,3,False)

# def print_primes(f):
#     print( [f[p] for p in primes[:20]] )

# print_primes(to_int(conv(p0,invert(p3))))
# print_primes(to_int(conv(p0,invert(p2))))
# print_primes(to_int(conv(p0,invert(p1))))
# print()
# print_primes(to_int(conv(p1,invert(p3))))
# print_primes(to_int(conv(p1,invert(p2))))
# print()
# print_primes(to_int(conv(p2,invert(p3))))


A = [0,2,6]
u = upsilon_gen(A,3,False)
r = rho_gen(A,3,False)
print(u[:40])
print(phi_gen(A)[:40])
print([u[i]-r[i] for i in range(40)])






