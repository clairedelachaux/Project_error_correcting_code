
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
    points = [(x[i], r[i]) for i in range(n)]
    return FX.lagrange_polynomial(points)

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
    if Q0%Q1 !=0:
        return "fail"
    f = -Q0//Q1
    return f.coefficients(sparse=False)


def  T(j, l, w, n, q):
    if l+j<w or l+w<j or w+j<l:
        return 0
    else:
        return sum(binomial(w, z)*binomial(w-z, z+j-l)*(q-2)**(l+w-j-2*z)*binomial(n-w, n-j-z)*(q-1)**(j+z-w) for z in range(0, w+1))

def A(w, n, k, q):
    if w == 0:
        return 1
    elif w < n-k+1:
        return 0
    else:
        return binomial(n, w)*sum((-1)**j*binomial(w, j)*(q**(w-(n-k)-j)-1) for j in range(0, w-(n-k+1)+1))

def P_error(n , k, q, p):
    t = (n-k)//2
    return sum((p/(q-1))**j*(1-p)**(n-j)*sum(A(w, n, k, q)*sum(T(j, l, w, n, q) for l in range (0, t+1)) for w in range (1, n+1)) for j in range (t+1, n+1))


def sim(iterations, k):
    # counters of errors and failures
    errors = 0
    failures = 0
    p = 0.007976
    P = [p, 1-p]
    error_dist=GeneralDiscreteDistribution(P)
    
    for i in range(iterations):
        x = vector(F, [a**(i-1) for i in range(1, n+1)])
        # create random message
        m = vector(F, [F.random_element() for i in range(k)])

        # encode message
        c = encode(m, x, k)

        # create error with computed error distribution and recieved word
        e = vector(F, [F.random_element() if error_dist.get_random_element()==0 else 0 for i in range(n)])
        if i % 500 == 0:
            print("iteration:", i, "\t decoding failures so far:", failures, "\t decoding errors so far:", errors)
        r = c + e

        # decode recieved word
        f = decode(r, x, n, k)

        # check if sent word and recieved word is the same
        if not list(m) == list(f)+[0]*(len(m)-len(list(f))):
            failures += 1
            if not f == "fail": 
                errors += 1
                print("Oh no! Error at iteration", i)
                print("\t Hammi+ng weight of current error vector: ", e.hamming_weight())
            else: 
                print("Phew! Just a failure at iteration", i)
                print("\t Hamming weight of current error vector: ", e.hamming_weight())
    return [failures, errors]

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
k = 245

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
n_sim = 10000        # Number of simulation
n_fail = 0       # Number of simulation where decoding failure occurs
t = (n - k)//2

# Perform simulations
l = sim(n_sim, k)
print(l)
