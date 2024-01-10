N = 2**129
L = 2*N**2

a = 2**128
while a % 24 != 1:
    a += 1

params1024 = []

from tqdm import tqdm
for _ in tqdm(range(10000000000)):
    a += 24
    if not is_pseudoprime(a):
        continue
    d = L - a^2
    assert d > 0
    #f = 4*(a^4 - 14*a^2*d + d^2)*(3*a^2 - d)*(a^2 - d)*(a^2 - 3*d)*a = 8*(a^4 - 14*a^2*d + d^2)*(3*a^2 - d)*(a^2 - 3*d)(a - N)(a + N)*a
    f1 = int(a - N)
    if f1 % 5 == 0:
        f1 //= 5
    f2 = int((a + N)/3)  
    if f2 % 5 == 0:
        f2 //= 5  
    f3 = int((3*a^2 - d)/4)
    if f3 % 5 == 0:
        f3 //= 5
    f4 = int((a^2 - 3*d)/4)  
    if f4 % 5 == 0:
        f4 //= 5 


    #print("-------")
    #print(factor(f1, limit=100)) 
    #print(factor(f2, limit=100)) 
    #print(factor(f3, limit=100)) 
    #print(factor(f4, limit=100)) 

    if all([is_pseudoprime(abs(f_fac)) for f_fac in [a, f1, f2, f3, f4]]):
        print("They were all prime")
        if all([kronecker_symbol(-d, abs(part)) == 1 for part in [f3, f4]]):

            K.<sq_d> = NumberField(x^2 + d)
            w = a + sq_d
            omega = w^12
            _, f = list(omega)

            print("Params:")
            print(f"d = {d}")
            print(f"d = {factor(d, limit = 10000)}")
            print(f"omega = {omega}")
            print(f"N(omega) = {factor(int(omega.norm()), limit=100)}")
            #assert f == 8*3*4*4*5*(a^4 - 14*a^2*d + d^2)*f1*f2*f3*f4*a
            #print(f"f = {f} = 8*3*4*4*{a^4 - 14*a^2*d + d^2}*{f1}*{f2}*{f3}*{f4}*{a}")
            print(f"f = {f} = junk*{f1}*{f2}*{f3}*{f4}*{a}")
            print(f"size d: {round(log(d, 2), 5)}")
            print(f"size a: {round(log(a, 2), 5)}")
            print(f"size f1: {round(log(abs(f1), 2), 5)}")
            print(f"size f2: {round(log(abs(f2), 2), 5)}")
            print(f"size f3: {round(log(abs(f3), 2), 5)}")
            print(f"size f4: {round(log(abs(f4), 2), 5)}")
            print("factoriszations:")
            print(f"d = {factor(d, limit = 10000)}")



