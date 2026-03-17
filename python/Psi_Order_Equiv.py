from base import *

mu = invert( [1]*bound )
La = conv(mu, eval(lambda n: math.log(n))) # Von Mangoldt
for i in range(bound):
    if abs(La[i]) < 1e-14: La[i] = 0.0

def standard_La_A_gen(A):
    def out(n):
        return math.prod(La[n+k] if n+k < bound else 1 for k in A)
    return eval( out )

A = [0,10]
singular = singular_series_of_A(A)
s_La = standard_La_A_gen( A )

Psi = sum(La) # Singular Psi
s_Psi_A = sum(s_La) # Standard Psi_A
x = bound

print( "Psi_A - C_A x", s_Psi_A - (x*singular))
print( "Psi_A - C_A Psi", s_Psi_A - (Psi*singular))