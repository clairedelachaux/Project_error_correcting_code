from main.sage import *

def T(j, l, w, n, q):
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


n=16
k=10
p=0.05
q=16
t = (n-k)//2

P1 = P_error(n, k, q, p)

# print(P1)

def sim(iterations, k):
    # counters of errors and failures
    errors = 0
    failures = 0

    # create code with k as input
    [encoder, decoder] = createCode(k)
    
    for i in range(iterations):
        # create random message
        m = vector(F, [F.random_element() for i in range(k)])

        # encode message
        c = encoder(m)

        # create error with computed error distribution and recieved word
        e = vector(F, [F.random_element() if error_dist.get_random_element()==0 else 0 for i in range(n)])
        if i % 500 == 0:
            print("iteration:", i, "\t decoding failures so far:", failures, "\t decoding errors so far:", errors)
        r = c + e

        # decode recieved word
        f = decoder(r)

        # check if sent word and recieved word is the same
        if not list(m) == list(f)+[0]*(len(m)-len(list(f))):
            failures += 1
            if not f == "fail": 
                errors += 1
                print("Oh no! Error at iteration", i)
                print("\t Hamming weight of current error vector: ", e.hamming_weight())
            else: 
                print("Phew! Just a failure at iteration", i)
                print("\t Hamming weight of current error vector: ", e.hamming_weight())
    return [failures, errors]


