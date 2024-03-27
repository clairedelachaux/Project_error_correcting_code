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
    return f.coefficients(sparse=False)




"""
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
"""


## Question 4
# Set global variables 
global F
global FX

# Set known parameters
q = 256
n = 255
k = 246

# Define a finite field with q = 256 elements and define a as an element of the finite field
a = chr(0x03B1)
F.<a> = FiniteField(q)

# Set a as the generator of the field
a = F.gen()

# Create polynomial ring FX with variable X where coefficients will be elements of F
FX.<X> = PolynomialRing(F)

# Generate the x vector of length n = 255 (covers all non-zero elements of our finite field using a)
x = vector(F, [a**(i-1) for i in range(1, n+1)])


# Perform simulations 
n_sim = 5        # Number of simulation
n_fail = 0       # Number of simulation where decoding failure occurs
t = (n - k)//2

# Perform simulations
for i in range(n_sim):
    # Generate a random message of length k
    m = vector([F.random_element() for _ in range(k)])

    # Encode the message
    r = encode(m,x,k)

    # Decode the received message
    m_hat = decode(r, x, n, k)

    # Check for decoding failure
    t_hat = sum(1 for mi, m_hat_i in zip(m, m_hat) if mi != m_hat_i)
    if t_hat > t:
        n_fail += 1

# Compute the probability of a decoding failure
P_fail = n_fail / n_sim
print("Probability of decoding failure:", P_fail)


