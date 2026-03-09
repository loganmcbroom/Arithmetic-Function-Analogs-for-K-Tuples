from base import *

def dm_conv(f,g,n):
    out = []
    for k in range(n):
        total = 0
        # Sum over a+b = k mod n
        for x in range(n):
            total += f(x) * g((k+n-x)%n)
        out.append(total)
    return out

def indicator(A):
    def out(n):
        return 1 if any(a == n for a in A) else 0
    return out

def chi(modulus):
    def out(n):
        return 1 if math.gcd(n,modulus) == 1 else 0
    return out

def totient_spectrum(A, n):
    return dm_conv(chi(n), indicator(A), n)

A = [0,1,3]
for n in range(1,13):
    print(totient_spectrum(A,n))
            

