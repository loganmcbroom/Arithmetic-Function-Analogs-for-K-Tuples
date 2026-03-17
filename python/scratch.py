from base import *

# primorial = 1
# A = [0]
# for p in primes[:10]:
#     primorial *= p
#     A.append(primorial)
#     print( p, 1.0/singular_series_of_A(A) )
print( singular_series_of_A([0]) )
print( singular_series_of_A([0,2]) )
print( singular_series_of_A([0,2,6]) )
print( singular_series_of_A([0,2,6,8]) )
print( singular_series_of_A([0,2,6,8,12]) )
print( singular_series_of_A([0,2,6,8,12,18]) )
print( singular_series_of_A([0,2,6,8,12,18,20]) )
print( singular_series_of_A([0,2,6,8,12,18,20,26]) )


# A = set([0,2,6])
# u_A = u_gen(A)
# w_A = w_gen(A)
# mu_A = mu_gen(A)
# phi_A = phi_gen(A)

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

# kts = [
#     [0],
#     [0,2],
#     [0,6],
#     [0,10],
#     [0,14],
#     [0,30],
#     [0,2,6],
#     [0,4,10],
#     [0,4,28],
#     [0,2,6,8],
#     [0,18,30,60],
#     [0,60,120,180],
#     [0,2,6,8,30],
#     ]

# uAs = [u_gen(A) for A in kts]
# wAs = [w_gen(A) for A in kts]

# index = 11
# 9, 10, 11 have prime 61??
# print( len() )
# print( conv_power(mu, 3)[:32] )
# u_cancel = conv( uAs[index], conv_power( mu, len(kts[index]) ) )
# w_cancel = conv(wAs[index], mu)

# print( sum(1 for x in u_cancel if x != 0) / bound )
# print( sum(1 for x in w_cancel if x != 0) / bound )

# Compute residue counts
# def r_gen(A): # u for primes
#     def out(p):
#         A_mod_p = set()
#         for k in A:
#             A_mod_p.add( k % p )
#         return len(A_mod_p)
#     return out
# Find modulated primes for each uA
# from primes import primes
# def find_mods(A, rA):
#     mods = []
#     for p in primes[:1000]:
#         if rA(p) < len(A):
#             mods.append(p)
#     return mods

# mods = find_mods( kts[index], r_gen(kts[index]) )
# print( mods )
# cancel_primes = []
# for p in primes[:100]:
#     for i in range( 1, bound ):
#         if u_cancel[i] != 0 and u_cancel[i] % p == 0: 
#             cancel_primes.append( p )
#             break
# print( cancel_primes )

# import pandas as pd
# print(pd.DataFrame({
#     "u" : u_cancel[:50],
#     # "w" : w_cancel[:50],
# }))


