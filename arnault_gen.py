from sage.all import *
from Crypto.Util.number import *
import random
def miller_rabin(n, b):
    """
    Miller Rabin test testing over all
    prime basis < b
    """
    basis = generate_basis(b)
    if n == 2 or n == 3:
        return True

    if n % 2 == 0:
        return False

    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for b in basis:
        x = pow(b, s, n)

        if x == 1 or x == n - 1:
            #print(b)
            continue
        #print(n-1)


        for _ in range(r - 1):
            x = pow(x, 2, n)
            #print("x:",x)
            if x == n - 1:
                #print("abc")
                break
        else:
            #print(b)
            return False
    return True

def generate_basis(n):
    basis = [True] * n
    for i in range(3, int(n**0.5)+1, 2):
        if basis[i]:
            basis[i*i::2*i] = [False]*((n-i*i-1)//(2*i)+1)
    return [2] + [i for i in range(3, n, 2) if basis[i]]

#
def find_pseudoprime(limit_primes,coeff,numbit_p1=5,search_range_bit=20):
    bases=generate_basis(limit_primes)
    #print(bases)
    list_a=[b * 4 for b in bases]
    #print(list_a)

    list_Sa=[]
    for a in bases:
        S_a=[]
        for i in range(1,a*4):
            #print("i:",i)
            if GCD(a*4,i)!=1:
                continue
            
            if kronecker(a,i)==-1:
                S_a.append(i)
        #print(f"base {a}, Sa: {S_a}")

        #print(f"base {a}")
        S_a_pr=set()
        for k in range(1,len(coeff)):
            S_a_i=set()
            if k==1:
                S_a_pr=set(S_a)
            
            for i in range(len(S_a)):
                S_a_i.add((S_a[i] + coeff[k] - 1) * inverse(coeff[k],4*a) % (4*a))

            S_a_pr = S_a_pr.intersection(S_a_i)
        #print(f"S_a_pr:{S_a_pr}")
        list_Sa.append(sorted(S_a_pr))
    #print(list_Sa)

    can_be_1 = all(any(x % 4 == 1 for x in opts) for opts in list_Sa)
    can_be_3 = all(any(x % 4 == 3 for x in opts) for opts in list_Sa)

    target_mod = None
    if can_be_3: 
        target_mod = 3
    elif can_be_1: 
        target_mod = 1
    else:
        print("Error: No consistent modulo 4 solution exists!")
        return
    
    #list_crt_a=[random.choice(x) for x in list_Sa]
    list_crt_a=[]
    for i in range(len(list_Sa)):
        while True:
            pick = random.choice(list_Sa[i]) #Randomly pick candidates
            if pick % 4 == target_mod:
                list_crt_a.append(pick)
                break
    #print(list_crt_a)
    #list_crt_a=[3, 7, 3, 15, 23, 47, 31, 47, 47, 55]
    list_crt_mod=list_a
    k2,k3=coeff[1],coeff[2]

    #print("k2:",k2)
    #print("k3:",k3)
    #print(inverse(-k3,k2),inverse(-k2,k3))
    list_crt_a.append(inverse(-k2,k3)%k3)
    list_crt_mod.append(k3)

    list_crt_a.append(inverse(-k3,k2)%k2)
    list_crt_mod.append(k2)

    try:
        base_p1=crt(list_crt_a,list_crt_mod)
        M=lcm(list_crt_mod)
        #print(base_p1)
        #print(M)
        
        base_p1_bit=int(base_p1).bit_length() # Base_p1 bits length 
        #Safety setting
        if numbit_p1<base_p1_bit:
            numbit_p1=base_p1_bit
        
        step=0 # step = pow(2,0) = 1
        if search_range_bit > 32:
            step += search_range_bit - 32 # Ensure always run in <= 2**32 iterations

        #Brute forcing p1 given p1 = z mod M
        for i in range(pow(2,numbit_p1-base_p1_bit),pow(2,numbit_p1 - base_p1_bit + search_range_bit),pow(2,step)):
            p1=i*M + base_p1
            if not isPrime(int(p1)):
                continue
            p2=k2*(p1-1) + 1
            if not isPrime(int(p2)):
                continue
            p3=k3*(p1-1) + 1
            if not isPrime(int(p3)):
                continue

            if miller_rabin(int(p1*p2*p3),limit_primes):
                print(f"p:",p1*p2*p3)
                print(f"p1:",p1)
                print(f"p2:",p2)
                print(f"p3:",p3)
                return (p1*p2*p3)#,[p1,p2,p3]
    except Exception as e:
        print("Error",e)




if __name__ == '__main__':
    """
    Pretty messy function for generating p1,p2,p3 whereas Carmichael number p1*p2*p3 = p can bypass 
    all prime bases of miller-rabin primality from 1 to n
    Ref: https://eprint.iacr.org/2018/749.pdf 


    numbit_p1: Minimum bit of p1
    limit_primes: Gen p that can bypass all prime from 1 up to number you choose 
    e.g limit_primes = 64 -> bypass [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
    coeff: coefficients of k1,k2 and k3
    k1 is always 1 while the rest need to be odd primes > limit_primes
    """
    #print('abc')
    numbit_p1=200
    search_range_bit=100
    limit_primes=64
    coeff=[1,101,313]
    p=find_pseudoprime(limit_primes,coeff,numbit_p1,search_range_bit)

    #print(p)
