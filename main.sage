def Van(x, k):
    n = len(x)
    M = matrix(F, k, n)
    for i in range(k):
        for j in range(n):
            M[i, j] = x[j]^i
    return M

def encode(m, x, k):
    V = Van(x, k)
    return m*V

def make_G(x):
    G = 1
    for xi in x:
        G *= (X-xi)
    return G

def L(x, i, n):
    xi = x[i]
    L = 1
    for h in range(i):
        L*=(X-x[h])/(x[i]-x[h])
    for h in range(i+1, n):
        L*=(X-x[h])/(x[i]-x[h])
    return L

def make_R(r, x, n):
    R=0
    for i in range(n):
        Li = L(x, i, n)
        R += r[i]*Li
    return R

def is_k_reduced(P, k):
    return P[0, 0].degree() >= P[0, 1].degree()+k-1 and P[1, 0].degree() < P[1, 1].degree()+k-1

def Q(P, k):
    degk_row_1 = max(P[0, 0].degree(), P[0, 1].degree()+k-1) 
    degk_row_2 = max(P[1, 0].degree(), P[1, 1].degree()+k-1)
    if degk_row_1<degk_row_2:
        return(P[0, 0], P[0, 1])
    return (P[1, 0], P[1, 1])

def decode(r, x, n, k):
    G = make_G(x)
    R = make_R(r, x, n)
    P = matrix(FX, [[G, 0], [-R, 1]])
    tau = floor(n-k/2)
    while not is_k_reduced(P, k):
        q = P[0, 0]//P[1, 0]
        r = P[0, 0]%P[1, 0]
        M = matrix(FX, [[0, 1], [1, -q]])
        P = M*P
    Q0, Q1 = Q(P, k)
    f = -Q0//Q1
    return f.coefficients()



##define the finite field
a = chr(0x03B1)
global F
global FX
q = 11
F.<a> = FiniteField(q)
a = F.gen()
FX.<X> = PolynomialRing(F)



n = 6
#x = vector(F, [a**(i-1) for i in range(1, n+1)])
x = vector(F, [1, 2, 3, 4, 5, 6])
k = 2
m = vector(F, [2, 1])
r = vector(F, [3, 4, 1, 6, 6, 8])
G = make_G(x)
R = make_R(r, x, n)
P = decode(r, x, n, k)
Chan = channels.StaticErrorRateChannel(F, 2)
print(Chan)

