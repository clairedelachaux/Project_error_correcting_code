

# This file was *autogenerated* from the file main.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_256 = Integer(256); _sage_const_255 = Integer(255); _sage_const_246 = Integer(246); _sage_const_0x03B1 = Integer(0x03B1); _sage_const_20 = Integer(20)
def Van(x, k):
    n = len(x)
    M = matrix(F, k, n)
    for i in range(k):
        for j in range(n):
            M[i, j] = x[j]**i
    return M

def encode(m, x, k):
    V = Van(x, k)
    return m*V

def make_G(x):
    G = _sage_const_1 
    for xi in x:
        G *= (X-xi)
    return G

def L(x, i, n):
    xi = x[i]
    L = _sage_const_1 
    for h in range(i):
        L*=(X-x[h])/(x[i]-x[h])
    for h in range(i+_sage_const_1 , n):
        L*=(X-x[h])/(x[i]-x[h])
    return L

def make_R(r, x, n):
    R=_sage_const_0 
    for i in range(n):
        Li = L(x, i, n)
        R += r[i]*Li
    return R

def is_k_reduced(P, k):
    return P[_sage_const_0 , _sage_const_0 ].degree() >= P[_sage_const_0 , _sage_const_1 ].degree()+k-_sage_const_1  and P[_sage_const_1 , _sage_const_0 ].degree() < P[_sage_const_1 , _sage_const_1 ].degree()+k-_sage_const_1 

def Q(P, k):
    degk_row_1 = max(P[_sage_const_0 , _sage_const_0 ].degree(), P[_sage_const_0 , _sage_const_1 ].degree()+k-_sage_const_1 ) 
    degk_row_2 = max(P[_sage_const_1 , _sage_const_0 ].degree(), P[_sage_const_1 , _sage_const_1 ].degree()+k-_sage_const_1 )
    if degk_row_1<degk_row_2:
        return(P[_sage_const_0 , _sage_const_0 ], P[_sage_const_0 , _sage_const_1 ])
    return (P[_sage_const_1 , _sage_const_0 ], P[_sage_const_1 , _sage_const_1 ])

def decode(r, x, n, k):
    G = make_G(x)
    R = make_R(r, x, n)
    P = matrix(FX, [[G, _sage_const_0 ], [-R, _sage_const_1 ]])
    tau = floor(n-k/_sage_const_2 )
    while not is_k_reduced(P, k):
        q = P[_sage_const_0 , _sage_const_0 ]//P[_sage_const_1 , _sage_const_0 ]
        r = P[_sage_const_0 , _sage_const_0 ]%P[_sage_const_1 , _sage_const_0 ]
        M = matrix(FX, [[_sage_const_0 , _sage_const_1 ], [_sage_const_1 , -q]])
        P = M*P
    Q0, Q1 = Q(P, k)
    f = -Q0//Q1
    return f.coefficients()





##define the finite field
#a = chr(0x03B1)
#global F
#global FX
#q = 11
#F.<a> = FiniteField(q)
#a = F.gen()
#FX.<X> = PolynomialRing(F)



#n = 6
#x = vector(F, [a**(i-1) for i in range(1, n+1)])
#x = vector(F, [1, 2, 3, 4, 5, 6])
#k = 2
#m = vector(F, [2, 1])
#r = vector(F, [3, 4, 1, 6, 6, 8])
#G = make_G(x)
#R = make_R(r, x, n)
#P = decode(r, x, n, k)
#Chan = channels.StaticErrorRateChannel(F, 2)
#print(Chan)


## Question 4
# Set global variables 
global F
global FX

# Set known parameters
q = _sage_const_256 
n = _sage_const_255 
k = _sage_const_246 

# Define a finite field with q = 256 elements and define a as an element of the finite filed
a = chr(_sage_const_0x03B1 )
F = FiniteField(q, names=('a',)); (a,) = F._first_ngens(1)

# Set a as the generator of the field
a = F.gen()

# Create polynomial ring FX with variable X where coefficients will be elements of F
FX = PolynomialRing(F, names=('X',)); (X,) = FX._first_ngens(1)

# Generate the x vector of length n = 255 (covers all non-zero elements of our finite field using a)
x = vector(F, [a**(i-_sage_const_1 ) for i in range(_sage_const_1 , n+_sage_const_1 )])


# Perform simulation 
n_sim = _sage_const_20        # Number of simulation
n_fail = _sage_const_0       # Number of simulation where decoding failure occurs
t = (n - k)//_sage_const_2 

for i in range(n_sim):
    # Generate a random message of length k
    m = vector([F.random_element() for _ in range(k)])
    r = encode(m,x,k)
    m_hat = decode(r, x, n, k)
    t_hat = sum(_sage_const_1  for mi, m_hat_i in zip(m, m_hat) if mi != m_hat_i)
    if t_hat > t:
        n_fail += _sage_const_1 

print(n_fail)



